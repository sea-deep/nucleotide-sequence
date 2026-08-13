import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { Alignment } from '../src/Alignment';

describe('Alignment module', () => {
  it('should perform Smith-Waterman local alignment', () => {
    // Expected:
    // query:     TGTTACGG
    // reference: GGTTGACTA
    // Local match should find GTT-AC or something similar depending on penalties
    // Let's use a simpler known case
    const q = new Seq().read('ACACACTA');
    const r = new Seq().read('AGCACACA');
    
    // With match=2, mismatch=-1, gapOpen=-2, gapExtend=-1
    // The best local alignment is 'ACACACTA' vs 'AGCACACA' ? 
    // Actually the 'ACACA' part matches perfectly.
    const res = Alignment.smithWaterman(q, r, { match: 2, mismatch: -1, gapOpen: -2, gapExtend: -1 });
    // q: ACACA
    // r: ACACA
    expect(res.queryAligned).toBe('ACACA');
    expect(res.referenceAligned).toBe('ACACA');
    expect(res.score).toBe(10); // 5 matches * 2 = 10
  });

  it('should handle gaps in Smith-Waterman', () => {
    const q = new Seq().read('ACACACTA');
    const r = new Seq().read('ACACTA');
    const res = Alignment.smithWaterman(q, r, { match: 2, mismatch: -1, gapOpen: -2, gapExtend: -1 });
    
    // Depending on gap penalty, it could be ACACACTA vs ACA--CTA
    // ACACTA matches exactly with ACACTA in q, score = 6 * 2 = 12
    expect(res.score).toBe(12);
  });

  it('should perform Needleman-Wunsch global alignment', () => {
    const q = new Seq().read('GATTACA');
    const r = new Seq().read('GCATGCG');
    // Global alignment should force end-to-end alignment
    const res = Alignment.needlemanWunsch(q, r, { match: 1, mismatch: -1, gapOpen: -1, gapExtend: -1 });
    
    // G-ATTACA
    // GCA-TGCG
    expect(res.queryAligned.length).toEqual(res.referenceAligned.length);
  });
});
