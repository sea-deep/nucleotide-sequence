import { Seq } from './Seq';
import { Alignment, AlignmentResult, AlignmentOptions } from './Alignment';

export interface HSP {
  queryStart: number;
  refStart: number;
  length: number;
  score: number;
}

/**
 * Implements a BLAST-like Seed-and-Extend heuristic search.
 * Much faster than full Smith-Waterman for large references.
 */
export class SeedSearch {
  private reference: Seq;
  private k: number;
  // 2-bit packed k-mer -> array of reference positions
  private index: Map<number, number[]>;

  constructor(reference: Seq, k: number = 11) {
    if (k > 15 || k < 3) throw new Error('k must be between 3 and 15 for 32-bit integer packing');
    this.reference = reference;
    this.k = k;
    this.index = new Map();
  }

  /**
   * Packs a k-mer into a 32-bit integer using 2-bit encoding.
   * A=00, C=01, G=10, T/U=11.
   * Returns -1 if the k-mer contains ambiguous bases (like N) which we don't index.
   */
  private pack(data: Uint8Array, start: number, k: number): number {
    let packed = 0;
    for (let i = 0; i < k; i++) {
      const byte = data[start + i];
      let val = 0;
      if (byte === 65) val = 0; // A
      else if (byte === 67) val = 1; // C
      else if (byte === 71) val = 2; // G
      else if (byte === 84 || byte === 85) val = 3; // T or U
      else return -1; // Ambiguous base, skip this seed

      packed = (packed << 2) | val;
    }
    return packed;
  }

  /**
   * Builds the k-mer hash index for the reference sequence.
   */
  public buildIndex(): this {
    this.index.clear();
    const data = this.reference['data'];
    for (let i = 0; i <= data.length - this.k; i++) {
      const packed = this.pack(data, i, this.k);
      if (packed !== -1) {
        const posArray = this.index.get(packed);
        if (posArray) {
          posArray.push(i);
        } else {
          this.index.set(packed, [i]);
        }
      }
    }
    return this;
  }

  /**
   * Searches the reference for the query sequence.
   * 1. Finds exact k-mer seeds.
   * 2. Performs ungapped extension to form HSPs.
   * 3. (Optional) Performs Smith-Waterman gapped alignment around HSPs.
   */
  public search(query: Seq, minScore: number = 22, gapped: boolean = true): AlignmentResult[] {
    const qData = query['data'];
    const rData = this.reference['data'];
    const hsps: HSP[] = [];
    const visitedDiagonals = new Set<number>();

    // 1 & 2. Seeding and Ungapped Extension
    for (let i = 0; i <= qData.length - this.k; i++) {
      const packed = this.pack(qData, i, this.k);
      if (packed === -1) continue;

      const hits = this.index.get(packed);
      if (hits) {
        for (const hit of hits) {
          // Diagonal = reference_pos - query_pos
          // Prevents extending the same HSP multiple times from overlapping seeds
          const diagonal = hit - i;
          if (visitedDiagonals.has(diagonal)) continue;
          visitedDiagonals.add(diagonal);

          // Ungapped extension left
          let left = 0;
          let score = this.k * 2; // Exact match seed = k * +2
          while (i - left - 1 >= 0 && hit - left - 1 >= 0 && qData[i - left - 1] === rData[hit - left - 1]) {
            left++;
            score += 2;
          }
          
          // Ungapped extension right
          let right = 0;
          while (
            i + this.k + right < qData.length && 
            hit + this.k + right < rData.length && 
            qData[i + this.k + right] === rData[hit + this.k + right]
          ) {
            right++;
            score += 2;
          }

          if (score >= minScore) {
            hsps.push({
              queryStart: i - left,
              refStart: hit - left,
              length: this.k + left + right,
              score
            });
          }
        }
      }
    }

    if (!gapped) {
      // Just return HSPs as pseud-alignments
      return hsps.sort((a, b) => b.score - a.score).map(hsp => {
        const qStr = String.fromCharCode.apply(null, Array.from(qData.subarray(hsp.queryStart, hsp.queryStart + hsp.length)));
        const rStr = String.fromCharCode.apply(null, Array.from(rData.subarray(hsp.refStart, hsp.refStart + hsp.length)));
        return {
          score: hsp.score,
          queryAligned: qStr,
          referenceAligned: rStr,
          queryStart: hsp.queryStart,
          referenceStart: hsp.refStart
        };
      });
    }

    // 3. Gapped Extension via Smith-Waterman
    const results: AlignmentResult[] = [];
    const PAD = Math.min(50, Math.floor(qData.length / 2)); 

    for (const hsp of hsps) {
      const qStart = Math.max(0, hsp.queryStart - PAD);
      const qEnd = Math.min(qData.length, hsp.queryStart + hsp.length + PAD);
      const rStart = Math.max(0, hsp.refStart - PAD);
      const rEnd = Math.min(rData.length, hsp.refStart + hsp.length + PAD);

      // Create temporary Seq objects for the local neighborhood
      const qSlice = new Seq();
      qSlice['data'] = qData.subarray(qStart, qEnd);
      
      const rSlice = new Seq();
      rSlice['data'] = rData.subarray(rStart, rEnd);

      const align = Alignment.smithWaterman(qSlice, rSlice);
      
      // Adjust start coordinates to global coordinates
      align.queryStart += qStart;
      align.referenceStart += rStart;
      
      results.push(align);
    }
    
    // Sort by score descending and return
    return results.sort((a, b) => b.score - a.score);
  }
}
