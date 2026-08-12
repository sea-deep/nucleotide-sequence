/**
 * The standard genetic code translation table.
 * Maps DNA codons (3 nucleotides) to Amino Acids.
 */
export const CODON_TABLE: Record<string, string> = {
  // Phenylalanine
  TTT: 'F', TTC: 'F',
  // Leucine
  TTA: 'L', TTG: 'L', CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  // Isoleucine
  ATT: 'I', ATC: 'I', ATA: 'I',
  // Methionine (Start)
  ATG: 'M',
  // Valine
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  // Serine
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S', AGT: 'S', AGC: 'S',
  // Proline
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  // Threonine
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  // Alanine
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  // Tyrosine
  TAT: 'Y', TAC: 'Y',
  // Stop Codons
  TAA: '*', TAG: '*', TGA: '*',
  // Histidine
  CAT: 'H', CAC: 'H',
  // Glutamine
  CAA: 'Q', CAG: 'Q',
  // Asparagine
  AAT: 'N', AAC: 'N',
  // Lysine
  AAA: 'K', AAG: 'K',
  // Aspartic Acid
  GAT: 'D', GAC: 'D',
  // Glutamic Acid
  GAA: 'E', GAG: 'E',
  // Cysteine
  TGT: 'C', TGC: 'C',
  // Tryptophan
  TGG: 'W',
  // Arginine
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R', AGA: 'R', AGG: 'R',
  // Glycine
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
};

import { Seq } from './Seq';

export class Translation {
  /**
   * Translates a sequence starting from a specific frame (0, 1, or 2).
   */
  public static translateFrame(seq: Seq, frame: number = 0): string {
    const sequenceStr = seq.sequence().toUpperCase();
    const normalizedSeq = seq.type === 'RNA' ? sequenceStr.replace(/U/g, 'T') : sequenceStr;
    
    let protein = '';
    for (let i = frame; i < normalizedSeq.length - 2; i += 3) {
      const codon = normalizedSeq.slice(i, i + 3);
      protein += CODON_TABLE[codon] || '?';
    }
    return protein;
  }

  /**
   * Translates a DNA/RNA sequence into an Amino Acid protein sequence.
   * Note: This strictly translates frame 0.
   */
  public static translate(seq: Seq): string {
    return this.translateFrame(seq, 0);
  }

  /**
   * Scans all 6 reading frames (3 forward, 3 reverse complement) for Open Reading Frames (ORFs).
   * Finds sequences that start with Methionine (M) and end with a Stop Codon (*).
   */
  public static findOpenReadingFrames(seq: Seq): string[] {
    const proteins: string[] = [];
    const orfRegex = /M[^*]*\*/g; // M followed by anything except *, ending with *
    
    const findInSeq = (s: Seq) => {
      for (let frame = 0; frame < 3; frame++) {
        const translated = this.translateFrame(s, frame);
        let match;
        while ((match = orfRegex.exec(translated)) !== null) {
          proteins.push(match[0].slice(0, -1)); // Remove the * stop codon from the output string
        }
      }
    };

    // Forward strand
    findInSeq(seq);
    
    // Reverse complement strand
    const compData = seq.complement().sequence();
    const revCompData = compData.split('').reverse().join('');
    const revSeq = new Seq(seq.type).read(revCompData);
    findInSeq(revSeq);
    
    return proteins;
  }
}
