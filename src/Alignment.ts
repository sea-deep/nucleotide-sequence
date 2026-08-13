import { Seq } from './Seq';

export interface AlignmentOptions {
  match?: number;
  mismatch?: number;
  gapOpen?: number;
  gapExtend?: number;
}

export interface AlignmentResult {
  score: number;
  queryAligned: string;
  referenceAligned: string;
  queryStart: number;
  referenceStart: number;
}

export class Alignment {
  /**
   * Smith-Waterman local alignment algorithm.
   * Finds the optimal local alignment between two sequences using dynamic programming
   * with affine gap penalties.
   */
  public static smithWaterman(
    query: Seq,
    reference: Seq,
    options: AlignmentOptions = {}
  ): AlignmentResult {
    const match = options.match ?? 2;
    const mismatch = options.mismatch ?? -1;
    const gapOpen = options.gapOpen ?? -2;
    const gapExtend = options.gapExtend ?? -1;

    const qStr = query.sequence().toUpperCase();
    const rStr = reference.sequence().toUpperCase();
    const m = qStr.length;
    const n = rStr.length;

    // H: main score matrix
    // E: gap in reference (insertion)
    // F: gap in query (deletion)
    const H = new Int32Array((m + 1) * (n + 1));
    const E = new Int32Array((m + 1) * (n + 1));
    const F = new Int32Array((m + 1) * (n + 1));

    let maxScore = 0;
    let maxI = 0;
    let maxJ = 0;

    for (let i = 1; i <= m; i++) {
      for (let j = 1; j <= n; j++) {
        const idx = i * (n + 1) + j;
        const upIdx = (i - 1) * (n + 1) + j;
        const leftIdx = i * (n + 1) + (j - 1);
        const diagIdx = (i - 1) * (n + 1) + (j - 1);

        E[idx] = Math.max(
          H[leftIdx] + gapOpen + gapExtend,
          E[leftIdx] + gapExtend
        );

        F[idx] = Math.max(
          H[upIdx] + gapOpen + gapExtend,
          F[upIdx] + gapExtend
        );

        const s = qStr[i - 1] === rStr[j - 1] ? match : mismatch;
        
        const score = Math.max(
          0,
          H[diagIdx] + s,
          E[idx],
          F[idx]
        );

        H[idx] = score;

        if (score > maxScore) {
          maxScore = score;
          maxI = i;
          maxJ = j;
        }
      }
    }

    // Backtrack
    let qAlign = '';
    let rAlign = '';
    let i = maxI;
    let j = maxJ;

    while (i > 0 && j > 0 && H[i * (n + 1) + j] > 0) {
      const idx = i * (n + 1) + j;
      const diagIdx = (i - 1) * (n + 1) + (j - 1);

      const s = qStr[i - 1] === rStr[j - 1] ? match : mismatch;

      if (H[idx] === H[diagIdx] + s) {
        qAlign = qStr[i - 1] + qAlign;
        rAlign = rStr[j - 1] + rAlign;
        i--;
        j--;
      } else if (H[idx] === E[idx]) {
        qAlign = '-' + qAlign;
        rAlign = rStr[j - 1] + rAlign;
        j--;
      } else if (H[idx] === F[idx]) {
        qAlign = qStr[i - 1] + qAlign;
        rAlign = '-' + rAlign;
        i--;
      } else {
        // Should not reach here for optimal path
        break;
      }
    }

    return {
      score: maxScore,
      queryAligned: qAlign,
      referenceAligned: rAlign,
      queryStart: i,
      referenceStart: j
    };
  }

  /**
   * Needleman-Wunsch global alignment algorithm.
   * Finds the optimal global alignment between two sequences using dynamic programming
   * with affine gap penalties.
   */
  public static needlemanWunsch(
    query: Seq,
    reference: Seq,
    options: AlignmentOptions = {}
  ): AlignmentResult {
    const match = options.match ?? 2;
    const mismatch = options.mismatch ?? -1;
    const gapOpen = options.gapOpen ?? -2;
    const gapExtend = options.gapExtend ?? -1;

    const qStr = query.sequence().toUpperCase();
    const rStr = reference.sequence().toUpperCase();
    const m = qStr.length;
    const n = rStr.length;

    const H = new Int32Array((m + 1) * (n + 1));
    const E = new Int32Array((m + 1) * (n + 1));
    const F = new Int32Array((m + 1) * (n + 1));

    // Initialize the boundaries
    const NEG_INF = -100000000;
    for (let i = 0; i <= m; i++) {
      H[i * (n + 1)] = i === 0 ? 0 : gapOpen + i * gapExtend;
      E[i * (n + 1)] = NEG_INF;
      F[i * (n + 1)] = i === 0 ? NEG_INF : gapOpen + i * gapExtend;
    }
    for (let j = 0; j <= n; j++) {
      H[j] = j === 0 ? 0 : gapOpen + j * gapExtend;
      E[j] = j === 0 ? NEG_INF : gapOpen + j * gapExtend;
      F[j] = NEG_INF;
    }

    for (let i = 1; i <= m; i++) {
      for (let j = 1; j <= n; j++) {
        const idx = i * (n + 1) + j;
        const upIdx = (i - 1) * (n + 1) + j;
        const leftIdx = i * (n + 1) + (j - 1);
        const diagIdx = (i - 1) * (n + 1) + (j - 1);

        E[idx] = Math.max(
          H[leftIdx] + gapOpen + gapExtend,
          E[leftIdx] + gapExtend
        );

        F[idx] = Math.max(
          H[upIdx] + gapOpen + gapExtend,
          F[upIdx] + gapExtend
        );

        const s = qStr[i - 1] === rStr[j - 1] ? match : mismatch;
        
        H[idx] = Math.max(
          H[diagIdx] + s,
          E[idx],
          F[idx]
        );
      }
    }

    // Backtrack
    let qAlign = '';
    let rAlign = '';
    let i = m;
    let j = n;

    while (i > 0 || j > 0) {
      const idx = i * (n + 1) + j;
      const diagIdx = i > 0 && j > 0 ? (i - 1) * (n + 1) + (j - 1) : 0;
      
      const s = (i > 0 && j > 0) ? (qStr[i - 1] === rStr[j - 1] ? match : mismatch) : 0;

      if (i > 0 && j > 0 && H[idx] === H[diagIdx] + s) {
        qAlign = qStr[i - 1] + qAlign;
        rAlign = rStr[j - 1] + rAlign;
        i--;
        j--;
      } else if (j > 0 && H[idx] === E[idx]) {
        qAlign = '-' + qAlign;
        rAlign = rStr[j - 1] + rAlign;
        j--;
      } else if (i > 0 && H[idx] === F[idx]) {
        qAlign = qStr[i - 1] + qAlign;
        rAlign = '-' + rAlign;
        i--;
      }
    }

    return {
      score: H[m * (n + 1) + n],
      queryAligned: qAlign,
      referenceAligned: rAlign,
      queryStart: 0,
      referenceStart: 0
    };
  }
}
