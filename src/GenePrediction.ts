import { Seq } from './Seq';

export interface GeneModel {
  start: number;
  end: number;
  strand: '+' | '-';
  frame: number;
  score: number;
  sequence: string;
}

/**
 * A simplified prokaryotic gene predictor inspired by Prodigal.
 * 1. Self-trains a codon usage model from long ORFs.
 * 2. Scores candidate ORFs using log-likelihood of codon bias.
 * 3. Applies a bonus for upstream Shine-Dalgarno (RBS) motifs.
 * 4. Greedily resolves overlaps to produce a final gene model.
 */
export class GenePrediction {
  private seq: Seq;
  private seqStr: string;
  private revSeqStr: string;

  private codingCodonFreq: Map<string, number> = new Map();
  private backgroundCodonFreq: Map<string, number> = new Map();

  constructor(seq: Seq) {
    this.seq = seq;
    this.seqStr = seq.sequence().toUpperCase();
    this.revSeqStr = seq.reverseComplement().sequence().toUpperCase();
  }

  /**
   * Main entry point to predict genes.
   */
  public predict(minCodons: number = 50): GeneModel[] {
    const candidates = this.findCandidateORFs(minCodons);
    if (candidates.length === 0) return [];

    this.trainModels(candidates);
    this.scoreCandidates(candidates);

    // Filter out very low scoring candidates
    const viable = candidates.filter(c => c.score > 0);

    return this.resolveOverlaps(viable);
  }

  /**
   * Finds all structural ORFs on both strands.
   */
  private findCandidateORFs(minCodons: number): GeneModel[] {
    const candidates: GeneModel[] = [];
    const minLength = minCodons * 3;
    const startCodons = ['ATG', 'GTG', 'TTG'];
    const stopCodons = ['TAA', 'TAG', 'TGA'];

    const searchStrand = (str: string, strand: '+' | '-') => {
      const len = str.length;
      for (let frame = 0; frame < 3; frame++) {
        let i = frame;
        while (i < len - 2) {
          const codon = str.substring(i, i + 3);
          if (startCodons.includes(codon)) {
            let j = i + 3;
            let foundStop = false;
            while (j < len - 2) {
              const nextCodon = str.substring(j, j + 3);
              if (stopCodons.includes(nextCodon)) {
                if (j - i >= minLength) {
                  candidates.push({
                    start: strand === '+' ? i : len - j - 3,
                    end: strand === '+' ? j + 3 : len - i,
                    strand,
                    frame: strand === '+' ? frame + 1 : -(frame + 1),
                    score: 0,
                    sequence: str.substring(i, j + 3)
                  });
                }
                i = j; // jump to stop codon
                foundStop = true;
                break;
              }
              j += 3;
            }
            if (!foundStop) break; // no more stops in this frame
          }
          i += 3;
        }
      }
    };

    searchStrand(this.seqStr, '+');
    searchStrand(this.revSeqStr, '-');

    return candidates;
  }

  /**
   * Self-trains codon usage. Assumes long ORFs (>300 codons) are highly likely to be true genes.
   */
  private trainModels(candidates: GeneModel[]): void {
    let totalBg = 0;
    let totalCoding = 0;

    // Background: just scan the whole sequence
    for (let i = 0; i < this.seqStr.length - 2; i += 3) {
      const codon = this.seqStr.substring(i, i + 3);
      this.backgroundCodonFreq.set(codon, (this.backgroundCodonFreq.get(codon) || 0) + 1);
      totalBg++;
    }

    // Coding: only from long ORFs
    const longORFs = candidates.filter(c => c.sequence.length >= 900);
    const trainingSet = longORFs.length > 0 ? longORFs : candidates; // fallback if no long ORFs

    for (const orf of trainingSet) {
      for (let i = 0; i < orf.sequence.length - 2; i += 3) {
        const codon = orf.sequence.substring(i, i + 3);
        this.codingCodonFreq.set(codon, (this.codingCodonFreq.get(codon) || 0) + 1);
        totalCoding++;
      }
    }

    // Normalize to probabilities
    for (const [codon, count] of this.backgroundCodonFreq.entries()) {
      this.backgroundCodonFreq.set(codon, count / totalBg);
    }
    for (const [codon, count] of this.codingCodonFreq.entries()) {
      this.codingCodonFreq.set(codon, count / totalCoding);
    }
  }

  /**
   * Scores each candidate based on Log-Likelihood Ratio of codon bias + RBS.
   */
  private scoreCandidates(candidates: GeneModel[]): void {
    for (const candidate of candidates) {
      let logLikelihood = 0;

      for (let i = 0; i < candidate.sequence.length - 2; i += 3) {
        const codon = candidate.sequence.substring(i, i + 3);
        const pCoding = Math.max(1e-5, this.codingCodonFreq.get(codon) || 1e-5);
        const pBg = Math.max(1e-5, this.backgroundCodonFreq.get(codon) || 1e-5);
        logLikelihood += Math.log(pCoding / pBg);
      }

      // Add RBS (Shine-Dalgarno) bonus
      const rbsScore = this.scoreRBS(candidate);
      
      // Weight the log-likelihood by length to avoid tiny ORFs dominating
      candidate.score = logLikelihood + rbsScore;
    }
  }

  /**
   * Looks for AGGAGG-like motifs 4-15nt upstream of the start codon.
   */
  private scoreRBS(candidate: GeneModel): number {
    let upstream = '';
    if (candidate.strand === '+') {
      const upStart = Math.max(0, candidate.start - 15);
      upstream = this.seqStr.substring(upStart, candidate.start - 4);
    } else {
      // For reverse strand, we need the 15nt BEFORE the start on the reverse strand
      // The start on the reverse strand is `candidate.end` in forward coords
      // But it's easier to just use revSeqStr
      const revStart = this.seqStr.length - candidate.end;
      const upStart = Math.max(0, revStart - 15);
      upstream = this.revSeqStr.substring(upStart, revStart - 4);
    }

    if (upstream.includes('AGGAGG')) return 15.0;
    if (upstream.includes('GGAGG') || upstream.includes('AGGAG')) return 10.0;
    if (upstream.includes('GGAG') || upstream.includes('GAGG')) return 5.0;
    
    return 0.0;
  }

  /**
   * Greedily resolves overlapping ORFs by taking the highest-scoring ones first.
   */
  private resolveOverlaps(candidates: GeneModel[]): GeneModel[] {
    // Sort by score descending
    candidates.sort((a, b) => b.score - a.score);
    
    const resolved: GeneModel[] = [];
    
    for (const candidate of candidates) {
      let overlap = false;
      for (const accepted of resolved) {
        // If they overlap on the same strand by a significant amount (e.g. > 30nt)
        if (candidate.strand === accepted.strand) {
          const overlapLength = Math.max(0, Math.min(candidate.end, accepted.end) - Math.max(candidate.start, accepted.start));
          if (overlapLength > 30) {
            overlap = true;
            break;
          }
        }
      }
      if (!overlap) {
        resolved.push(candidate);
      }
    }

    // Return sorted by genomic position
    return resolved.sort((a, b) => a.start - b.start);
  }
}
