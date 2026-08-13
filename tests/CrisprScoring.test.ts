import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { CrisprScoring } from '../src/CrisprScoring';

describe('CrisprScoring', () => {
  it('should find spacers and calculate on-target scores', () => {
    // Spacer: GGCGAGGAGCTGTTCACCGG (20) + NGG
    const seq = new Seq().read('GGCGAGGAGCTGTTCACCGGTGG');
    const targets = CrisprScoring.findSpacers(seq, 'NGG');
    
    expect(targets.length).toBe(1);
    expect(targets[0].spacer).toBe('GGCGAGGAGCTGTTCACCGG');
    expect(targets[0].onTargetScore).toBeDefined();
    // G at position 20 gives a boost
    expect(targets[0].onTargetScore).toBeGreaterThan(50);
  });
  
  it('should penalize poly-T runs in on-target score', () => {
    const score = CrisprScoring.calculateOnTargetScore('GGCGAGGAGCTTTTTACCGG');
    expect(score).toBeLessThan(50);
  });
  
  it('should calculate CFD scores correctly', () => {
    const guide = 'GGCGAGGAGCTGTTCACCGG';
    
    // Perfect match
    expect(CrisprScoring.calculateCFDScore(guide, guide)).toBe(1.0);
    
    // Mismatch in seed region (pos 18) -> heavy penalty
    const seedMismatch = 'GGCGAGGAGCTGTTCACCGT'; // G->T at pos 20
    expect(CrisprScoring.calculateCFDScore(guide, seedMismatch)).toBeCloseTo(0.1);
    
    // Mismatch in PAM-distal region (pos 3) -> light penalty
    const distalMismatch = 'GGAGAGGAGCTGTTCACCGG'; // C->A at pos 3
    expect(CrisprScoring.calculateCFDScore(guide, distalMismatch)).toBeCloseTo(0.9);
  });
});
