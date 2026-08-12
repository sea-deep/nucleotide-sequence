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

  it('should compute GC content correctly', () => {
    const seq = new Seq().read('GCGCATAT');
    expect(seq.gcContent()).toBe(0.5); // 4 GC out of 8
  });

  it('should generate reverse complement correctly', () => {
    const seq = new Seq().read('ATGC');
    const rc = seq.reverseComplement();
    expect(rc.sequence()).toBe('GCAT');
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
    // W(A/T)->W, S(C/G)->S, M(A/C)->K(T/G), K->M, R(A/G)->Y(T/C), Y->R, B->V, V->B, D->H, H->D, N->N
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
    expect(content.GC).toBe(0.5); // 50% GC
  });

  it('should calculate melting temperature (Tm) correctly', () => {
    // Primer: ATGC ATGC ATGC ATGC (8 A/T, 8 G/C)
    // Tm = 2*(8) + 4*(8) = 16 + 32 = 48
    const seq = new Seq().read('ATGCATGCATGCATGC');
    expect(seq.meltingTemperature()).toBe(48);
  });

  it('should correctly parse FASTA format', () => {
    const fasta = `>Sequence 1\nATGC\n>Comment\nCGTA`;
    const seq = new Seq().readFASTA(fasta);
    expect(seq.sequence()).toBe('ATGCCGTA');
  });

  it('should correctly calculate molecular weight', () => {
    const seq = new Seq().read('A');
    expect(seq.molecularWeight()).toBe(313.2);
  });

  it('should splice sequences correctly', () => {
    const seq = new Seq().read('ATGC');
    // Remove 2 elements from index 1 (TG) -> AC
    const spliced = seq.splice(1, 2);
    expect(spliced.sequence()).toBe('AC');
    
    // Insert 'TT' at index 1 -> ATTTGC
    const insertSeq = new Seq().read('TT');
    const inserted = seq.splice(1, 0, insertSeq);
    expect(inserted.sequence()).toBe('ATTTGC');
  });

  it('should correctly parse FASTQ format', () => {
    const fastq = `@SEQ_ID\nGATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65`;
    const seq = new Seq().readFASTQ(fastq);
    expect(seq.sequence()).toBe('GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT');
  });

  it('should correctly extract K-mers', () => {
    const seq = new Seq().read('ATGC');
    const k3 = seq.kmers(3);
    expect(k3).toEqual(['ATG', 'TGC']);
    const k2 = seq.kmers(2);
    expect(k2).toEqual(['AT', 'TG', 'GC']);
  });

  it('should find CRISPR targets (PAM sequences)', () => {
    // NGG pam -> any + GG
    const seq = new Seq().read('ATGG C TGG C AGG C CGG C');
    const targets = seq.findCRISPRTargets('NGG');
    // ATGG (PAM is TGG at 1), CTGG (PAM is TGG at 5), CAGG (PAM is AGG at 9), CCGG (PAM is CGG at 13)
    expect(targets).toEqual([1, 5, 9, 13]);
  });
});
