import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';

describe('Seq Class', () => {
  it('should read sequence data correctly', () => {
    const seq = new Seq();
    seq.read('ATGC');
    expect(seq.sequence()).toBe('ATGC');
    expect(seq.size()).toBe(4);
  });

  it('should ignore spaces and newlines', () => {
    const seq = new Seq();
    seq.read('AT GC\nTA');
    expect(seq.sequence()).toBe('ATGCTA');
  });

  it('should reject invalid characters', () => {
    expect(() => new Seq().read('ATG123')).toThrow(TypeError);
    expect(() => new Seq().read('ATGZ')).toThrow(TypeError);
    expect(() => new Seq().read('ATG!')).toThrow(TypeError);
  });

  it('should accept all IUPAC degenerate codes', () => {
    // Should not throw
    const seq = new Seq().read('ACGTUMRWSYKVHDBN');
    expect(seq.size()).toBe(16);
  });

  it('should compute GC content correctly', () => {
    const seq = new Seq().read('GCGCATAT');
    expect(seq.gcContent()).toBe(0.5);
  });

  it('should generate reverse complement correctly for standard bases', () => {
    const seq = new Seq().read('ATGC');
    const rc = seq.reverseComplement();
    expect(rc.sequence()).toBe('GCAT');
  });

  it('should generate reverse complement correctly for IUPAC degenerate codes', () => {
    // This was a bug: reverseComplement() previously only handled A/T/C/G.
    // M(A/C)->K(T/G), R(A/G)->Y(T/C) — complement then reverse.
    const seq = new Seq().read('ATGMR');
    const rc = seq.reverseComplement();
    // Complement: TACKM -> reversed: MKCAT
    // Wait, let me think carefully:
    // A->T, T->A, G->C, M->K, R->Y
    // Complement of ATGMR = TACKY
    // Reverse of TACKY = YKCAT
    expect(rc.sequence()).toBe('YKCAT');
  });

  it('should compute complementary DNA sequence correctly', () => {
    const seq = new Seq('DNA').read('ATGC');
    const comp = seq.complement();
    expect(comp.sequence()).toBe('TACG');
    expect(comp.type).toBe('DNA');
  });

  it('should compute complementary RNA sequence correctly', () => {
    const seq = new Seq('RNA').read('AUGC');
    const comp = seq.complement();
    expect(comp.sequence()).toBe('UACG');
    expect(comp.type).toBe('RNA');
  });

  it('should properly complement degenerate nucleotides', () => {
    const seq = new Seq().read('WSMKRYBVDHN');
    const comp = seq.complement();
    expect(comp.sequence()).toBe('WSKMYRVBHDN');
  });

  it('should transcribe DNA to RNA', () => {
    const dna = new Seq('DNA').read('ATGC');
    const rna = dna.transcribe();
    expect(rna.type).toBe('RNA');
    expect(rna.sequence()).toBe('AUGC');
  });

  it('should reverse transcribe RNA to DNA', () => {
    const rna = new Seq('RNA').read('AUGC');
    const dna = rna.reverseTranscribe();
    expect(dna.type).toBe('DNA');
    expect(dna.sequence()).toBe('ATGC');
  });

  it('should calculate GC content and fractional content correctly', () => {
    const seq = new Seq().read('ATGC');
    const content = seq.fractionalContent();
    expect(content.A).toBe(0.25);
    expect(content.G).toBe(0.25);
    expect(content.GC).toBe(0.5);
  });

  it('should calculate melting temperature (Tm) for valid-length primers', () => {
    // 16nt primer: ATGC repeated 4x (8 A/T, 8 G/C)
    // Tm = 2*(8) + 4*(8) = 48
    const seq = new Seq().read('ATGCATGCATGCATGC');
    expect(seq.meltingTemperature()).toBe(48);
  });

  it('should throw RangeError for Tm on sequences outside 14-20nt', () => {
    const shortSeq = new Seq().read('ATGCATGC'); // 8nt
    expect(() => shortSeq.meltingTemperature()).toThrow(RangeError);

    const longSeq = new Seq().read('ATGCATGCATGCATGCATGCATGCATGC'); // 28nt
    expect(() => longSeq.meltingTemperature()).toThrow(RangeError);
  });

  it('should correctly parse FASTA format', () => {
    const fasta = `>Sequence 1\nATGC\n>Comment\nCGTA`;
    const records = Seq.readFASTA(fasta);
    expect(records.length).toBe(2);
    expect(records[0].id).toBe('Sequence');
    expect(records[0].description).toBe('1');
    expect(records[0].seq.sequence()).toBe('ATGC');
    expect(records[1].id).toBe('Comment');
    expect(records[1].seq.sequence()).toBe('CGTA');
  });

  it('should correctly calculate molecular weight', () => {
    const seq = new Seq().read('A');
    expect(seq.molecularWeight()).toBeCloseTo(251.25, 2);
  });

  it('should splice sequences correctly', () => {
    const seq = new Seq().read('ATGC');
    const spliced = seq.splice(1, 2);
    expect(spliced.sequence()).toBe('AC');
    
    const insertSeq = new Seq().read('TT');
    const inserted = seq.splice(1, 0, insertSeq);
    expect(inserted.sequence()).toBe('ATTTGC');
  });

  it('should correctly parse FASTQ format', () => {
    const fastq = `@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>CCCCCCCC65`;
    const records = Seq.readFASTQ(fastq);
    expect(records.length).toBe(1);
    expect(records[0].id).toBe('SEQ_ID');
    expect(records[0].seq.sequence()).toBe('GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT');
    
    const decoder = new TextDecoder();
    expect(decoder.decode(records[0].quality)).toBe(`!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>CCCCCCCC65`);
  });

  it('should correctly extract K-mers', () => {
    const seq = new Seq().read('ATGC');
    const decoder = new TextDecoder();
    
    const k3 = Array.from(seq.kmers(3)).map(arr => decoder.decode(arr));
    expect(k3).toEqual(['ATG', 'TGC']);
    
    const k2 = Array.from(seq.kmers(2)).map(arr => decoder.decode(arr));
    expect(k2).toEqual(['AT', 'TG', 'GC']);
  });

  it('should find PAM sites', () => {
    const seq = new Seq().read('ATGGCTGGCAGGCCGG');
    const sites = seq.findPAMSites('NGG');
    // Should find positions where ?GG occurs
    expect(sites.length).toBeGreaterThan(0);
    // Verify each found site actually matches ?GG
    const seqStr = seq.sequence();
    for (const pos of sites) {
      expect(seqStr[pos + 1]).toBe('G');
      expect(seqStr[pos + 2]).toBe('G');
    }
  });

  it('should calculate precise molecular weight using standard IDT constants', () => {
    // ATCG residues: A(313.21) + T(304.20) + C(289.18) + G(329.21) = 1235.8
    const seq = new Seq().read('ATCG');
    
    // Standard linear synthetic oligo (5'-OH, 3'-OH) 
    // Expect: 1235.8 - 61.96 = 1173.84
    const stdWeight = seq.molecularWeight();
    expect(stdWeight).toBeCloseTo(1173.84, 2);
    
    // Phosphorylated oligo (5'-PO4, 3'-OH)
    // Expect: 1235.8 + 18.02 = 1253.82
    const phosWeight = seq.molecularWeight({ phosphorylated: true });
    expect(phosWeight).toBeCloseTo(1253.82, 2);
  });
});
