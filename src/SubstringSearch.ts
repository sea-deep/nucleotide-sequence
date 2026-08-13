import { Seq, SeqType } from './Seq';

export class MatchResult {
  constructor(
    public readonly position: number,
    public readonly matches: number,
    private readonly sequenceStr: string,
    public readonly type: SeqType
  ) {}

  /**
   * Returns a Seq object representing the matched portion of the reference sequence.
   */
  public matchedRegion(): Seq {
    return new Seq(this.type).read(this.sequenceStr);
  }
}

/**
 * SubstringSearch provides ungapped substring matching with wildcard tolerance.
 *
 * For each position in the reference, it counts the number of matching
 * nucleotides against the query (Hamming distance complement). 'N' is
 * treated as a universal wildcard on both query and reference.
 *
 * This is NOT sequence alignment. It does not use gap penalties,
 * substitution matrices, or dynamic programming. For gapped alignment,
 * use Smith-Waterman or Needleman-Wunsch.
 */
export class SubstringSearch {
  private query: Seq;
  private reference: Seq;
  private matchScores: Uint32Array;

  constructor(query: Seq, reference: Seq) {
    if (query.type !== reference.type) {
      throw new Error('Query and Reference sequences must be of the same type (DNA or RNA)');
    }
    this.query = query;
    this.reference = reference;
    const possiblePositions = Math.max(0, reference.size() - query.size() + 1);
    this.matchScores = new Uint32Array(possiblePositions);
  }

  /**
   * Execute the search. Computes match scores at every position.
   */
  public initialize(): this {
    // We access the private 'data' via bracket notation to avoid changing Seq's public API
    const qData = this.query['data'] as Uint8Array;
    const rData = this.reference['data'] as Uint8Array;
    const qLen = qData.length;
    const rLen = rData.length;
    const N_CODE = 78; // 'N'

    for (let i = 0; i <= rLen - qLen; i++) {
      let score = 0;
      for (let j = 0; j < qLen; j++) {
        const qChar = qData[j];
        const rChar = rData[i + j];
        if (qChar === rChar || qChar === N_CODE || rChar === N_CODE) {
          score++;
        }
      }
      this.matchScores[i] = score;
    }
    return this;
  }

  /**
   * Returns the top results sorted by match count.
   * Uses a bounded selection (min-heap) for O(N) performance
   * instead of sorting the entire array.
   */
  public top(maxResults: number = 5): MatchResult[] {
    const n = this.matchScores.length;
    if (n === 0) return [];

    // Collect top-k indices using a simple bounded insertion sort.
    // For small k (typical usage), this is faster than a heap.
    const topIndices: number[] = [];
    const topScores: number[] = [];

    for (let i = 0; i < n; i++) {
      const score = this.matchScores[i];
      if (score === 0) continue;

      if (topIndices.length < maxResults) {
        // Still filling the buffer — insert in sorted position
        let pos = topIndices.length;
        while (pos > 0 && topScores[pos - 1] < score) {
          topScores[pos] = topScores[pos - 1];
          topIndices[pos] = topIndices[pos - 1];
          pos--;
        }
        topScores[pos] = score;
        topIndices[pos] = i;
      } else if (score > topScores[topIndices.length - 1]) {
        // Replace the smallest entry
        let pos = topIndices.length - 1;
        while (pos > 0 && topScores[pos - 1] < score) {
          topScores[pos] = topScores[pos - 1];
          topIndices[pos] = topIndices[pos - 1];
          pos--;
        }
        topScores[pos] = score;
        topIndices[pos] = i;
      }
    }

    const rStr = this.reference.sequence();
    const qLen = this.query.size();
    return topIndices.map((idx, i) =>
      new MatchResult(idx, topScores[i], rStr.slice(idx, idx + qLen), this.reference.type)
    );
  }

  /**
   * Returns the single best MatchResult.
   */
  public best(): MatchResult | null {
    let bestIdx = -1;
    let maxScore = -1;
    for (let i = 0; i < this.matchScores.length; i++) {
      if (this.matchScores[i] > maxScore) {
        maxScore = this.matchScores[i];
        bestIdx = i;
      }
    }
    if (bestIdx === -1) return null;
    return new MatchResult(bestIdx, maxScore, this.reference.sequence().slice(bestIdx, bestIdx + this.query.size()), this.reference.type);
  }
}
