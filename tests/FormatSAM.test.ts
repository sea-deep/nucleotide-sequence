import { describe, it, expect } from 'vitest';
import { FormatSAM } from '../src/FormatSAM';

describe('FormatSAM', () => {
  it('should parse a standard SAM file', () => {
    const samData = `
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
read1	16	ref	7	30	8M2I4M1D3M	*	0	0	TTAGATAAAGGATACTG	*	NM:i:1	MD:Z:8^A7
`.trim();

    const records = FormatSAM.parse(samData);
    expect(records.length).toBe(1);
    
    const r = records[0];
    expect(r.qname).toBe('read1');
    expect(r.flag).toBe(16);
    expect(r.rname).toBe('ref');
    expect(r.pos).toBe(7);
    expect(r.mapq).toBe(30);
    expect(r.cigar).toBe('8M2I4M1D3M');
    expect(r.seq).toBe('TTAGATAAAGGATACTG');
    
    // Check tags
    expect(r.tags['NM']).toBe(1);
    expect(r.tags['MD']).toBe('8^A7');
  });
  
  it('should throw on invalid format', () => {
    const invalidData = `read1\t16\tref`;
    expect(() => FormatSAM.parse(invalidData)).toThrow(/expected at least 11 fields/);
  });
});
