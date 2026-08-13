import { describe, it, expect, beforeEach } from 'vitest';
import { Seq } from '../src/Seq';
import { Parallel } from '../src/Parallel';
import { Worker as NodeWorker } from 'worker_threads';

// Mock Web Worker environment for Node using real worker_threads!
class MockWorker {
  onmessage: ((e: any) => void) | null = null;
  onerror: ((e: any) => void) | null = null;
  private worker: NodeWorker;

  constructor(script: string) {
    const wrapper = `
      const { parentPort } = require('worker_threads');
      global.self = global;
      self.postMessage = (data, transfer) => parentPort.postMessage(data, transfer);
      parentPort.on('message', async (e) => {
        if (self.onmessage) {
          try {
             await self.onmessage({ data: e });
          } catch(err) {
             console.error(err);
          }
        }
      });
      ${script}
    `;
    this.worker = new NodeWorker(wrapper, { eval: true });
    this.worker.on('message', (msg) => {
      if (this.onmessage) this.onmessage({ data: msg });
    });
    this.worker.on('error', (err) => {
      if (this.onerror) this.onerror(err);
    });
  }
  postMessage(message: any) {
    this.worker.postMessage(message);
  }
  terminate() {
    this.worker.terminate();
  }
}

describe('Parallel module', () => {
  beforeEach(() => {
    // Inject Mocks for zeroworker
    global.Worker = MockWorker as any;
    global.Blob = class Blob {
      content: string;
      constructor(public parts: any[]) {
        this.content = parts.join('');
      }
    } as any;
    global.URL = { createObjectURL: (blob: any) => blob.content } as any;
  });

  it('should align sequences in parallel', async () => {
    const query = new Seq().read('ACACACTA');
    const ref1 = new Seq().read('AGCACACA');
    const ref2 = new Seq().read('ACACTA');
    
    const results = await Parallel.align(query, [ref1, ref2], { match: 2, mismatch: -1, gapOpen: -2, gapExtend: -1 });
    
    expect(results.length).toBe(2);
    expect(results[0].score).toBe(10); 
    expect(results[1].score).toBe(12);
  });

  it('should count kmers in parallel', async () => {
    const seq = new Seq().read('ATGCATGCATGC');
    // ATGC (4), k=4
    const counts = await Parallel.kmerCount(seq, 4, 2); // 2 chunks
    
    expect(counts['ATGC']).toBe(3);
    expect(counts['TGCA']).toBe(2);
    expect(counts['GCAT']).toBe(2);
    expect(counts['CATG']).toBe(2);
  });
});
