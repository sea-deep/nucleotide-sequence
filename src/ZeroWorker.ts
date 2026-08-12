import { Seq } from './Seq';
import { MatchResult } from './MatchMap';

const workerCode = `
self.onmessage = function(e) {
  const { type, query, reference, pam, k } = e.data;
  
  if (type === 'align') {
    let bestPosition = -1;
    let maxMatches = -1;
    for (let i = 0; i <= reference.length - query.length; i++) {
      let matches = 0;
      for (let j = 0; j < query.length; j++) {
        const q = query[j];
        const r = reference[i + j];
        if (q === r || q === 78 || r === 78) matches++;
      }
      if (matches > maxMatches) {
        maxMatches = matches;
        bestPosition = i;
      }
      if (matches === query.length) break;
    }
    self.postMessage({ type: 'align', position: bestPosition, matches: maxMatches });
  } 
  
  else if (type === 'crispr') {
    const targets = [];
    const pamLen = pam.length;
    for (let i = 0; i <= reference.length - pamLen; i++) {
      let match = true;
      for (let p = 0; p < pamLen; p++) {
        const pCode = pam.charCodeAt(p);
        const rCode = reference[i + p];
        if (pCode === 78) continue; // N
        if (pCode === 82 && (rCode === 65 || rCode === 71)) continue; // R
        if (pCode !== rCode) { match = false; break; }
      }
      if (match) targets.push(i);
    }
    self.postMessage({ type: 'crispr', targets });
  }
  
  else if (type === 'kmers') {
    const kMers = [];
    for (let i = 0; i <= reference.length - k; i++) {
      const slice = new Uint8Array(k);
      for (let j = 0; j < k; j++) slice[j] = reference[i + j];
      kMers.push(String.fromCharCode.apply(null, slice));
    }
    self.postMessage({ type: 'kmers', kmers: kMers });
  }
};
`;

/**
 * A highly optimized, zero-dependency Web Worker proxy for genomic algorithms.
 */
export class ZeroWorker {
  private worker: Worker | null = null;
  private isSupported: boolean = false;
  private messageId = 0;

  constructor() {
    if (typeof Worker !== 'undefined' && typeof Blob !== 'undefined' && typeof URL !== 'undefined') {
      const blob = new Blob([workerCode], { type: 'application/javascript' });
      this.worker = new Worker(URL.createObjectURL(blob));
      this.isSupported = true;
    }
  }

  private postAndWait(message: any, transferables?: Transferable[]): Promise<any> {
    return new Promise((resolve, reject) => {
      this.worker!.onmessage = (e) => resolve(e.data);
      this.worker!.onerror = (e) => reject(e);
      if (transferables) {
        this.worker!.postMessage(message, transferables);
      } else {
        this.worker!.postMessage(message);
      }
    });
  }

  public async align(querySeq: Seq, refSeq: Seq, transferOwnership: boolean = false): Promise<MatchResult> {
    if (!this.isSupported || !this.worker) {
      const { MatchMap } = await import('./MatchMap');
      const map = new MatchMap(querySeq, refSeq).initialize();
      return map.best()!;
    }

    const qArray = querySeq.raw();
    const rArray = refSeq.raw();
    const isShared = typeof SharedArrayBuffer !== 'undefined' && (qArray.buffer instanceof SharedArrayBuffer || rArray.buffer instanceof SharedArrayBuffer);

    let data;
    if (isShared) {
      data = await this.postAndWait({ type: 'align', query: qArray, reference: rArray });
    } else if (transferOwnership) {
      data = await this.postAndWait({ type: 'align', query: qArray, reference: rArray }, [qArray.buffer, rArray.buffer]);
    } else {
      data = await this.postAndWait({ type: 'align', query: qArray, reference: rArray });
    }

    return {
      position: data.position,
      matches: data.matches,
      alignment: () => {
        try {
          const alignedArray = rArray.slice(data.position, data.position + querySeq.length());
          const alignedSeq = new Seq(refSeq.type());
          (alignedSeq as any)['data'] = alignedArray;
          return alignedSeq;
        } catch (err) {
          throw new Error('Memory detached.');
        }
      }
    };
  }

  public async findCRISPRTargets(refSeq: Seq, pam: string = 'NGG'): Promise<number[]> {
    if (!this.isSupported || !this.worker) return refSeq.findCRISPRTargets(pam);
    const rArray = refSeq.raw();
    const data = await this.postAndWait({ type: 'crispr', reference: rArray, pam });
    return data.targets;
  }

  public async kmers(refSeq: Seq, k: number): Promise<string[]> {
    if (!this.isSupported || !this.worker) return refSeq.kmers(k);
    const rArray = refSeq.raw();
    const data = await this.postAndWait({ type: 'kmers', reference: rArray, k });
    return data.kmers;
  }

  public terminate() {
    if (this.worker) {
      this.worker.terminate();
      this.worker = null;
    }
  }
}
