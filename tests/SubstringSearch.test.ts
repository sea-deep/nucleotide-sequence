import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { SubstringSearch } from '../src/SubstringSearch';

describe('SubstringSearch', () => {
  it('should find the exact match index correctly', () => {
    const query = new Seq().read('ATGC');
    const reference = new Seq().read('TTTTTTATGCTTTTTT');
    
    const search = new SubstringSearch(query, reference).initialize();
    const best = search.best();
    expect(best?.position).toBe(6);
    expect(best?.matchedRegion().sequence()).toBe('ATGC');
  });

  it('should return top results sorted by matches', () => {
    const query = new Seq().read('ATGC');
    const reference = new Seq().read('ATTCTTTTATGC');
    
    const search = new SubstringSearch(query, reference).initialize();
    const top = search.top(2);
    
    expect(top.length).toBe(2);
    expect(top[0].position).toBe(8);
    expect(top[0].matches).toBe(4);
    expect(top[0].matchedRegion().sequence()).toBe('ATGC');
    
    expect(top[1].position).toBe(0);
    expect(top[1].matches).toBe(3); 
    expect(top[1].matchedRegion().sequence()).toBe('ATTC');
  });

  it('should support N wildcard matches', () => {
    const query = new Seq().read('ATNC');
    const reference = new Seq().read('TTTATGC');
    
    const search = new SubstringSearch(query, reference).initialize();
    const top = search.top(1);
    
    expect(top[0].position).toBe(3);
    expect(top[0].matches).toBe(4);
  });

  it('should return empty results for no matches', () => {
    const query = new Seq().read('AAAA');
    const reference = new Seq().read('CCCC');
    
    const search = new SubstringSearch(query, reference).initialize();
    const top = search.top(5);
    expect(top.length).toBe(0);
  });

  it('should throw if query and reference types differ', () => {
    const query = new Seq('DNA').read('ATGC');
    const reference = new Seq('RNA').read('AUGC');
    expect(() => new SubstringSearch(query, reference)).toThrow();
  });
});
