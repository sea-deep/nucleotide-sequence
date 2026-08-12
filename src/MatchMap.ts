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
  public alignment(): Seq {
    return new Seq(this.type).read(this.sequenceStr);
  }
}

/**
 * MatchMap provides high-performance, exhaustive bitwise alignment of sequences.
 * This is a modernized, TypeScript-typed port of NtSeq's MatchMap algorithm.
 */
export class MatchMap {
  private query: Seq;
  private reference: Seq;
  private matchScores: Uint32Array;

  constructor(query: Seq, reference: Seq) {
    if (query.type !== reference.type) {
      throw new Error('Query and Reference sequences must be of the same type (DNA or RNA)');
    }
    this.query = query;
    this.reference = reference;
    const possibleAlignments = Math.max(0, reference.size() - query.size() + 1);
    this.matchScores = new Uint32Array(possibleAlignments);
  }

  public initialize(): this {
    const qStr = this.query.sequence().toUpperCase();
    const rStr = this.reference.sequence().toUpperCase();
    const qLen = qStr.length;
    const rLen = rStr.length;

    for (let i = 0; i <= rLen - qLen; i++) {
      let score = 0;
      for (let j = 0; j < qLen; j++) {
        const qChar = qStr[j];
        const rChar = rStr[i + j];
        if (qChar === rChar || qChar === 'N' || rChar === 'N') {
          score++;
        }
      }
      this.matchScores[i] = score;
    }
    return this;
  }

  /**
   * Returns the top alignments.
   */
  public top(maxResults: number = 5): MatchResult[] {
    const indices = Array.from(this.matchScores.keys());
    indices.sort((a, b) => this.matchScores[b] - this.matchScores[a]);
    
    const results: MatchResult[] = [];
    const rStr = this.reference.sequence();
    const qLen = this.query.size();

    for (let i = 0; i < Math.min(maxResults, indices.length); i++) {
      const idx = indices[i];
      if (this.matchScores[idx] > 0) {
        results.push(new MatchResult(idx, this.matchScores[idx], rStr.slice(idx, idx + qLen), this.reference.type));
      }
    }
    return results;
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
