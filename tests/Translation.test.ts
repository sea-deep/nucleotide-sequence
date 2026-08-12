import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { Translation } from '../src/Translation';

describe('Translation Module', () => {
  it('should translate DNA to Amino Acids correctly', () => {
    // ATG (Methionine/Start) -> GCC (Alanine) -> TAA (Stop)
    const seq = new Seq().read('ATGCCCTAA');
    const protein = Translation.translate(seq);
    expect(protein).toBe('MP*'); // CCC is Proline, so ATG CCC TAA -> M P *
  });

  it('should translate RNA to Amino Acids correctly', () => {
    const seq = new Seq('RNA').read('AUGCCUUAA'); 
    // AUG -> M, CCU -> P, UAA -> *
    const protein = Translation.translate(seq);
    expect(protein).toBe('MP*');
  });

  it('should pad unknown codons with ?', () => {
    const seq = new Seq().read('ATGNNNTAA');
    const protein = Translation.translate(seq);
    expect(protein).toBe('M?*');
  });

  it('should find Open Reading Frames (ORFs)', () => {
    // ATGCCTTAA has ATG (M), CCT (P), TAA (*)
    // Let's add some junk before it: CCC ATGCCTTAA GGG
    const seq = new Seq().read('CCCATGCCTTAAGGG');
    const orfs = Translation.findOpenReadingFrames(seq);
    // Should find MP (Methionine, Proline)
    expect(orfs).toContain('MP');
  });
});
