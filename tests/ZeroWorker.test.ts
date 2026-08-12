import { describe, it, expect, afterEach } from 'vitest';
import { Seq } from '../src/Seq';
import { ZeroWorker } from '../src/ZeroWorker';

describe('ZeroWorker Multithreading Engine', () => {
  let worker: ZeroWorker;

  afterEach(() => {
    if (worker) {
      worker.terminate();
    }
  });

  it('should fallback to synchronous MatchMap if Worker is unavailable in environment', async () => {
    // In Vitest (Node.js), window.Worker is typically undefined unless polyfilled.
    // Our wrapper should gracefully fallback to the synchronous MatchMap!
    worker = new ZeroWorker();
    
    const ref = new Seq().read('ATGCCTGGATGC');
    const query = new Seq().read('CCTG');

    const result = await worker.align(query, ref);
    expect(result.position).toBe(3);
    expect(result.matches).toBe(4);
  });

  it('should correctly extract alignment using the fallback alignment() helper', async () => {
    worker = new ZeroWorker();
    
    const ref = new Seq().read('ATGCCTGGATGC');
    const query = new Seq().read('CCTG');

    const result = await worker.align(query, ref);
    expect(result.alignment().sequence()).toBe('CCTG');
  });

  it('should support degenerate wildcard matching in fallback mode', async () => {
    worker = new ZeroWorker();
    
    const ref = new Seq().read('ATGNCTGGATGC');
    const query = new Seq().read('CCTG'); // N matches C

    const result = await worker.align(query, ref);
    expect(result.position).toBe(3);
    expect(result.matches).toBe(4);
  });

  it('should compute kmers via worker fallback', async () => {
    worker = new ZeroWorker();
    const ref = new Seq().read('ATGC');
    const kmers = await worker.kmers(ref, 2);
    expect(kmers).toEqual(['AT', 'TG', 'GC']);
  });

  it('should compute CRISPR targets via worker fallback', async () => {
    worker = new ZeroWorker();
    const ref = new Seq().read('ATGCTGG');
    const targets = await worker.findCRISPRTargets(ref, 'TGG');
    expect(targets).toEqual([4]);
  });
});
