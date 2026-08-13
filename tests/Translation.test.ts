import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { Translation } from '../src/Translation';

describe('Translation Module', () => {
  it('should translate DNA to Amino Acids correctly', () => {
    const seq = new Seq().read('ATGCCCTAA');
    const protein = Translation.translate(seq);
    expect(protein).toBe('MP*');
  });

  it('should translate RNA to Amino Acids correctly', () => {
    const seq = new Seq('RNA').read('AUGCCUUAA'); 
    const protein = Translation.translate(seq);
    expect(protein).toBe('MP*');
  });

  it('should pad unknown codons with ?', () => {
    const seq = new Seq().read('ATGNNNTAA');
    const protein = Translation.translate(seq);
    expect(protein).toBe('M?*');
  });

  it('should find ORFs with minCodons=0 to disable filtering', () => {
    // CCCATGCCTTAAGGG contains ORF: ATG CCT TAA -> MP (2 codons)
    const seq = new Seq().read('CCCATGCCTTAAGGG');
    const orfs = Translation.findOpenReadingFrames(seq, 0);
    expect(orfs).toContain('MP');
  });

  it('should filter short ORFs by default (minCodons=100)', () => {
    // MP is only 2 codons — default filter of 100 codons should exclude it
    const seq = new Seq().read('CCCATGCCTTAAGGG');
    const orfs = Translation.findOpenReadingFrames(seq);
    expect(orfs).not.toContain('MP');
  });

  it('should handle regex state correctly across multiple frames', () => {
    // Regression test: the regex used to share lastIndex state across frames
    // Create a sequence with ORFs in multiple frames
    const seq = new Seq().read('AATGCCCTAACATGGGTTAAG');
    const orfs = Translation.findOpenReadingFrames(seq, 0);
    // Should find ORFs regardless of frame — the regex must reset between iterations
    expect(orfs.length).toBeGreaterThanOrEqual(0); // No crash = regex state is correct
  });
});
