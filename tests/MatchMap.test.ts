import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { MatchMap } from '../src/MatchMap';

describe('MatchMap Algorithm', () => {
  it('should find the exact match index correctly', () => {
    const query = new Seq().read('ATGC');
    const reference = new Seq().read('TTTTTTATGCTTTTTT');
    
    const map = new MatchMap(query, reference).initialize();
    const best = map.best();
    expect(best?.position).toBe(6);
    expect(best?.alignment().sequence()).toBe('ATGC');
  });

  it('should return top results sorted by matches', () => {
    const query = new Seq().read('ATGC');
    const reference = new Seq().read('ATTCTTTTATGC');
    
    const map = new MatchMap(query, reference).initialize();
    const top = map.top(2);
    
    expect(top.length).toBe(2);
    expect(top[0].position).toBe(8);
    expect(top[0].matches).toBe(4);
    expect(top[0].alignment().sequence()).toBe('ATGC');
    
    expect(top[1].position).toBe(0);
    expect(top[1].matches).toBe(3); 
    expect(top[1].alignment().sequence()).toBe('ATTC');
  });

  it('should support N wildcard matches', () => {
    const query = new Seq().read('ATNC'); // N can be anything
    const reference = new Seq().read('TTTATGC'); // G will match N
    
    const map = new MatchMap(query, reference).initialize();
    const top = map.top(1);
    
    expect(top[0].position).toBe(3);
    expect(top[0].matches).toBe(4); // A=A, T=T, N=G(wildcard), C=C
  });
});
