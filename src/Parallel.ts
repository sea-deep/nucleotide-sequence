import { Seq } from './Seq';
import { AlignmentOptions, AlignmentResult } from './Alignment';
import { WorkerPool } from 'zeroworker';

export class Parallel {
  private static alignmentPool: WorkerPool | null = null;
  private static kmerPool: WorkerPool | null = null;

  /**
   * Aligns a query sequence against an array of reference sequences in parallel.
   */
  public static async align(query: Seq, references: Seq[], options: AlignmentOptions = {}): Promise<AlignmentResult[]> {
    if (!this.alignmentPool) {
      // Inline Smith-Waterman logic for the worker since it cannot access outside scope.
      // This is necessary because ZeroWorker serializes the function via toString().
      const workerLogic = (payload: { qData: Uint8Array, rData: Uint8Array, match: number, mismatch: number, gapOpen: number, gapExtend: number }) => {
        const { qData, rData, match, mismatch, gapOpen, gapExtend } = payload;
        const m = qData.length;
        const n = rData.length;

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

            E[idx] = Math.max(H[leftIdx] + gapOpen + gapExtend, E[leftIdx] + gapExtend);
            F[idx] = Math.max(H[upIdx] + gapOpen + gapExtend, F[upIdx] + gapExtend);

            const s = qData[i - 1] === rData[j - 1] ? match : mismatch;
            
            const score = Math.max(0, H[diagIdx] + s, E[idx], F[idx]);
            H[idx] = score;

            if (score > maxScore) {
              maxScore = score;
              maxI = i;
              maxJ = j;
            }
          }
        }

        let qAlign = '';
        let rAlign = '';
        let i = maxI;
        let j = maxJ;

        while (i > 0 && j > 0 && H[i * (n + 1) + j] > 0) {
          const idx = i * (n + 1) + j;
          const diagIdx = (i - 1) * (n + 1) + (j - 1);
          const upIdx = (i - 1) * (n + 1) + j;
          const leftIdx = i * (n + 1) + (j - 1);

          const s = qData[i - 1] === rData[j - 1] ? match : mismatch;

          if (H[idx] === H[diagIdx] + s) {
            qAlign = String.fromCharCode(qData[i - 1]) + qAlign;
            rAlign = String.fromCharCode(rData[j - 1]) + rAlign;
            i--;
            j--;
          } else if (H[idx] === E[idx]) {
            qAlign = '-' + qAlign;
            rAlign = String.fromCharCode(rData[j - 1]) + rAlign;
            j--;
          } else if (H[idx] === F[idx]) {
            qAlign = String.fromCharCode(qData[i - 1]) + qAlign;
            rAlign = '-' + rAlign;
            i--;
          } else {
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
      };

      this.alignmentPool = new WorkerPool(workerLogic);
      await this.alignmentPool.initialize();
    }

    const qData = query['data'];
    const match = options.match ?? 2;
    const mismatch = options.mismatch ?? -1;
    const gapOpen = options.gapOpen ?? -2;
    const gapExtend = options.gapExtend ?? -1;

    const promises = references.map(ref => {
      // Create a copy of qData for each worker if we want to retain the original,
      // but for Zero-Copy, we should transfer. Since we need to reuse qData, we MUST clone it.
      // However, rData belongs to the reference sequence. If we transfer it, the reference Seq object
      // will lose its buffer. Therefore, we should only transfer a cloned buffer.
      // Zero-copy is only truly beneficial if we don't need the buffer on the main thread anymore.
      // To properly utilize ZeroWorker's speed without breaking Seq objects, we slice (copy) here,
      // but explicitly transfer the slice to the worker to avoid a SECOND structured clone penalty.
      const qDataClone = new Uint8Array(qData);
      const rDataClone = new Uint8Array(ref['data']);

      return this.alignmentPool!.execute({
        qData: qDataClone,
        rData: rDataClone,
        match,
        mismatch,
        gapOpen,
        gapExtend
      }, [qDataClone.buffer, rDataClone.buffer]);
    });

    return Promise.all(promises);
  }

  /**
   * Parallel k-mer counting.
   * Chunks a sequence into segments, counts k-mers in parallel, and merges the results.
   */
  public static async kmerCount(seq: Seq, k: number, chunks: number = 4): Promise<Record<string, number>> {
    if (k <= 0 || k > seq.size()) {
      throw new Error('Invalid k-mer size');
    }

    if (!this.kmerPool) {
      const workerLogic = (payload: { data: Uint8Array, k: number }) => {
        const { data, k } = payload;
        const counts: Record<string, number> = {};
        const decoder = new TextDecoder();
        
        for (let i = 0; i <= data.length - k; i++) {
          // Decode is slightly slower but works inside worker
          const kmer = decoder.decode(data.subarray(i, i + k));
          counts[kmer] = (counts[kmer] || 0) + 1;
        }
        return counts;
      };

      this.kmerPool = new WorkerPool(workerLogic);
      await this.kmerPool.initialize();
    }

    const data = seq['data'];
    const len = data.length;
    const segmentSize = Math.ceil(len / chunks);
    const promises: Promise<Record<string, number>>[] = [];

    for (let i = 0; i < chunks; i++) {
      const start = i * segmentSize;
      let end = start + segmentSize;
      
      // Need to overlap by k-1 to not miss k-mers across chunk boundaries
      if (i < chunks - 1) {
        end += k - 1;
      }
      
      if (start >= len) break;
      end = Math.min(end, len);

      const chunkData = new Uint8Array(data.subarray(start, end));

      promises.push(this.kmerPool.execute({
        data: chunkData,
        k
      }, [chunkData.buffer]));
    }

    const results = await Promise.all(promises);
    
    // Merge counts
    const finalCounts: Record<string, number> = {};
    for (const res of results) {
      for (const [kmer, count] of Object.entries(res)) {
        finalCounts[kmer] = (finalCounts[kmer] || 0) + count;
      }
    }
    
    return finalCounts;
  }
}
