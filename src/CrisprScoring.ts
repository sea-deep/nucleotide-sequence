import { Seq } from './Seq';

export interface CrisprTarget {
  position: number;
  strand: '+' | '-';
  spacer: string;       // 20nt sequence
  pam: string;          // PAM sequence
  gcContent: number;    // GC% of spacer (30-70% is ideal)
  onTargetScore?: number; // Doench '16 / Rule Set 2 (0-100)
}

export interface OffTargetHit {
  position: number;
  strand: '+' | '-';
  sequence: string;
  mismatches: number;
  cfdScore: number; // 0.0 to 1.0
}

/**
 * Static CFD weights for mismatches.
 * Rows: Position 1 to 20.
 * Cols: Mismatch type (e.g., rA:dC, rA:dG, etc.)
 * Simplified representative matrix for SpCas9 based on Doench et al.
 */
const CFD_WEIGHTS: Record<number, Record<string, number>> = {
  // Pos 1-20 (5' to 3')
  1:  { 'rA:dC': 1.0, 'rA:dG': 1.0, 'rA:dT': 1.0, 'rC:dA': 1.0, 'rC:dG': 1.0, 'rC:dT': 1.0, 'rG:dA': 1.0, 'rG:dC': 1.0, 'rG:dT': 1.0, 'rU:dA': 1.0, 'rU:dC': 1.0, 'rU:dG': 1.0 },
  // ... omitting full 20x12 matrix for brevity, using positional penalties
  // Positions closer to PAM (18-20) are penalized heavily (score ~0.1)
  // Positions far from PAM (1-5) are penalized lightly (score ~0.9)
};

// Simplified positional weight matrix for On-Target scoring (Rule Set 2 proxy)
// Weights for A, C, G, T at positions 1-20
const ON_TARGET_WEIGHTS: Record<number, Record<string, number>> = {
  20: { A: -0.1, C: -0.1, G: 0.5, T: -0.2 }, // G at pos 20 is highly favored
  16: { A: 0.2, C: 0.1, G: -0.1, T: -0.2 },
  // ... simplified proxy
};

export class CrisprScoring {
  /**
   * Finds all CRISPR target sites (spacers) in a sequence given a PAM.
   */
  public static findSpacers(seq: Seq, pam: string = 'NGG', spacerLength: number = 20): CrisprTarget[] {
    const targets: CrisprTarget[] = [];
    const seqStr = seq.sequence().toUpperCase();
    
    // Find on forward strand (+)
    const fwdPamSites = seq.findPAMSites(pam);
    for (const pamPos of fwdPamSites) {
      if (pamPos >= spacerLength) {
        const spacer = seqStr.substring(pamPos - spacerLength, pamPos);
        const gcCount = (spacer.match(/[GC]/g) || []).length;
        
        targets.push({
          position: pamPos - spacerLength,
          strand: '+',
          spacer,
          pam: seqStr.substring(pamPos, pamPos + pam.length),
          gcContent: gcCount / spacerLength,
          onTargetScore: this.calculateOnTargetScore(spacer)
        });
      }
    }
    
    // Find on reverse strand (-)
    const revSeq = seq.reverseComplement();
    const revSeqStr = revSeq.sequence().toUpperCase();
    const revPamSites = revSeq.findPAMSites(pam);
    const len = seq.size();
    
    for (const revPamPos of revPamSites) {
      if (revPamPos >= spacerLength) {
        const spacer = revSeqStr.substring(revPamPos - spacerLength, revPamPos);
        const gcCount = (spacer.match(/[GC]/g) || []).length;
        
        targets.push({
          position: len - revPamPos,
          strand: '-',
          spacer,
          pam: revSeqStr.substring(revPamPos, revPamPos + pam.length),
          gcContent: gcCount / spacerLength,
          onTargetScore: this.calculateOnTargetScore(spacer)
        });
      }
    }
    
    return targets;
  }

  /**
   * Calculates a simplified on-target efficiency score (0-100) based on positional nucleotide frequencies.
   */
  public static calculateOnTargetScore(spacer: string): number {
    if (spacer.length !== 20) return 0;
    
    // Start with a baseline score of 50
    let score = 50.0;
    
    // Penalize poly-T runs (terminates U6 transcription)
    if (spacer.includes('TTTT')) score -= 40;
    
    // GC content bounds (30% to 80% is ideal)
    const gcCount = (spacer.match(/[GC]/g) || []).length;
    const gcPercent = gcCount / 20;
    if (gcPercent < 0.3 || gcPercent > 0.8) score -= 15;
    
    // Apply positional weights
    for (let i = 0; i < 20; i++) {
      const pos = i + 1;
      const nuc = spacer[i];
      if (ON_TARGET_WEIGHTS[pos] && ON_TARGET_WEIGHTS[pos][nuc]) {
        score += ON_TARGET_WEIGHTS[pos][nuc] * 10;
      }
    }
    
    // G at position 20 is highly favored for U6 promoter efficiency
    if (spacer[19] === 'G') score += 10;
    if (spacer[19] === 'C') score -= 5;
    
    return Math.max(0, Math.min(100, Math.round(score * 10) / 10));
  }

  /**
   * Calculates the CFD off-target score between a guide and an off-target sequence.
   * Returns a value between 0.0 and 1.0 (1.0 = perfect match).
   */
  public static calculateCFDScore(guide: string, offTarget: string): number {
    if (guide.length !== offTarget.length || guide.length !== 20) return 0.0;
    
    let score = 1.0;
    for (let i = 0; i < 20; i++) {
      if (guide[i] !== offTarget[i]) {
        const pos = i + 1;
        // Simplified CFD penalty: mismatches in the "seed" region (pos 10-20) are penalized heavily
        // Mismatches in the PAM-distal region (pos 1-9) are penalized lightly
        let penalty = 1.0;
        if (pos >= 15) penalty = 0.1;
        else if (pos >= 10) penalty = 0.3;
        else if (pos >= 5) penalty = 0.7;
        else penalty = 0.9;
        
        score *= penalty;
      }
    }
    
    return score;
  }
}
