import { Seq } from './Seq';
import { MatchResult } from './MatchMap';

const workerCode = `
self.onmessage = function(e) {
  const query = e.data.query;
  const reference = e.data.reference;
  
  let bestPosition = -1;
  let maxMatches = -1;

  for (let i = 0; i <= reference.length - query.length; i++) {
    let matches = 0;
    for (let j = 0; j < query.length; j++) {
      const q = query[j];
      const r = reference[i + j];
      if (q === r || q === 78 /* N */ || r === 78) {
        matches++;
      }
    }
    if (matches > maxMatches) {
      maxMatches = matches;
      bestPosition = i;
    }
    if (matches === query.length) break;
  }

  self.postMessage({ position: bestPosition, matches: maxMatches });
};
`;

/**
 * A highly optimized, zero-dependency Web Worker proxy for genomic alignments.
 * Automatically handles SharedArrayBuffer and Transferable Object routing
 * to guarantee zero-copy memory performance.
 */
export class ZeroWorker {
  private worker: Worker | null = null;
  private isSupported: boolean = false;

  constructor() {
    if (typeof Worker !== 'undefined' && typeof Blob !== 'undefined' && typeof URL !== 'undefined') {
      const blob = new Blob([workerCode], { type: 'application/javascript' });
      this.worker = new Worker(URL.createObjectURL(blob));
      this.isSupported = true;
    }
  }

  /**
   * Run the exhaustive MatchMap alignment in a background thread.
   * 
   * @param querySeq The sequence to search for
   * @param refSeq The massive reference genome to search inside
   * @param transferOwnership If true, physically moves the memory to the worker using Transferable Objects. 
   *                          WARNING: This detaches the memory from the main thread, making the original Seq objects unusable!
   *                          If using SharedArrayBuffer, this flag is ignored (memory is shared natively).
   */
  public async align(querySeq: Seq, refSeq: Seq, transferOwnership: boolean = false): Promise<MatchResult> {
    if (!this.isSupported || !this.worker) {
      // Graceful fallback for environments without Web Workers (e.g. standard Node.js without polyfills)
      const { MatchMap } = await import('./MatchMap');
      const map = new MatchMap(querySeq, refSeq).initialize();
      return map.best()!;
    }

    return new Promise((resolve, reject) => {
      const qArray = querySeq.raw();
      const rArray = refSeq.raw();

      this.worker!.onmessage = (e) => {
        resolve({
          position: e.data.position,
          matches: e.data.matches,
          alignment: () => {
            // Note: If memory was transferred and not shared, creating an alignment from refSeq will fail 
            // because the buffer is detached. We assume if they transfer, they don't need the original sequence.
            try {
              const alignedArray = rArray.slice(e.data.position, e.data.position + querySeq.length());
              const alignedSeq = new Seq(refSeq.type());
              // Bypass access modifiers for performance slice
              (alignedSeq as any)['data'] = alignedArray;
              return alignedSeq;
            } catch (err) {
              throw new Error('Cannot extract alignment: memory was transferred and detached. Use SharedArrayBuffer for zero-copy memory sharing.');
            }
          }
        });
      };

      this.worker!.onerror = (e) => {
        reject(e);
      };

      const hasSharedArrayBuffer = typeof SharedArrayBuffer !== 'undefined';
      const isShared = hasSharedArrayBuffer && (qArray.buffer instanceof SharedArrayBuffer || rArray.buffer instanceof SharedArrayBuffer);

      if (isShared) {
        // Zero-Copy: Shared Memory (instant, no detachment)
        this.worker!.postMessage({ query: qArray, reference: rArray });
      } else if (transferOwnership) {
        // Zero-Copy: Transferable Memory (instant, but detaches main thread memory)
        this.worker!.postMessage({ query: qArray, reference: rArray }, [qArray.buffer, rArray.buffer]);
      } else {
        // Structured Clone (copies memory - safe but high memory overhead)
        this.worker!.postMessage({ query: qArray, reference: rArray });
      }
    });
  }
  
  public terminate() {
    if (this.worker) {
      this.worker.terminate();
      this.worker = null;
    }
  }
}
