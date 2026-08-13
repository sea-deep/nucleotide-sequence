import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { SeedSearch } from '../src/Search';

describe('SeedSearch', () => {
  it('should find seeds and perform gapped extension', () => {
    // Reference has a 21-nt region
    const ref = new Seq().read('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT');
    // Query has a match in the middle but with a small insertion "AAA"
    const query = new Seq().read('CGTACGTACGAAAACGTACGT');

    const searcher = new SeedSearch(ref, 5).buildIndex();
    const results = searcher.search(query, 10, true);

    expect(results.length).toBeGreaterThan(0);
    // Should match around score 2 * 18 - 2 = 34 ish
    expect(results[0].score).toBeGreaterThan(20);
  });
  
  it('should return ungapped HSPs when gapped=false', () => {
    const ref = new Seq().read('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT');
    const query = new Seq().read('CGTACGTACGAAAACGTACGT');

    const searcher = new SeedSearch(ref, 5).buildIndex();
    const results = searcher.search(query, 10, false);

    expect(results.length).toBeGreaterThan(0);
    // Should find the two exact pieces (length 10 and length 8)
    expect(results.some(r => r.score === 20)).toBe(true); // CGTACGTACG
    expect(results.some(r => r.score === 16)).toBe(true); // ACGTACGT
  });
});
