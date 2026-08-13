import { Seq } from './Seq';

export const TRANSLATION_TABLES: Record<number, Record<string, string>> = {
  // Table 1: Standard
  1: {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
  },
  // Table 2: Vertebrate Mitochondrial
  2: {
    'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': '*', 'AGG': '*',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W',
  },
  // Table 11: Bacterial, Archaeal and Plant Plastid
  11: {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
  }
};

export class Translation {
  /**
   * Translates a nucleotide sequence into an amino acid sequence.
   * By default uses NCBI Translation Table 1 (Standard).
   * Note: RNA sequences are internally converted (U->T) for codon lookup.
   */
  public static translate(seq: Seq, tableId: number = 1): string {
    return this.translateFrame(seq, 0, tableId);
  }

  /**
   * Translates a sequence starting from a specific frame offset (0, 1, or 2).
   */
  public static translateFrame(seq: Seq, frame: number = 0, tableId: number = 1): string {
    const table = TRANSLATION_TABLES[tableId];
    if (!table) {
      throw new Error(`Unsupported translation table ID: ${tableId}`);
    }

    const seqStr = seq.sequence().toUpperCase();
    
    let protein = '';
    for (let i = frame; i < seqStr.length - 2; i += 3) {
      const codon = seqStr.substring(i, i + 3).replace(/U/g, 'T');
      protein += table[codon] || '?'; // Use ? for ambiguous/degenerate codons
    }
    return protein;
  }


  /**
   * Scans all 6 reading frames to find structural Open Reading Frames (ORFs).
   * 
   * Note: This strictly finds start-to-stop sequences. It does NOT predict
   * biological translation viability (e.g., it ignores Ribosome Binding Sites 
   * or Kozak consensus sequences), and scanning raw eukaryotic genomic DNA 
   * will yield massive false positives due to introns.
   *
   * @param seq The input sequence to scan.
   * @param minCodons Minimum ORF length in codons (default: 100). Set to 0 to disable filtering.
   * @param tableId Translation table to use (default: 1).
   */
  public static findOpenReadingFrames(seq: Seq, minCodons: number = 100, tableId: number = 1): string[] {
    const proteins: string[] = [];
    
    const findInSeq = (s: Seq) => {
      for (let frame = 0; frame < 3; frame++) {
        const translated = this.translateFrame(s, frame, tableId);
        // Create a fresh regex per iteration to avoid stale lastIndex state
        const orfRegex = /M[^*]*\*/g;
        let match;
        while ((match = orfRegex.exec(translated)) !== null) {
          const orf = match[0].slice(0, -1); // Remove the * stop codon
          if (orf.length >= minCodons) {
            proteins.push(orf);
          }
        }
      }
    };

    // Forward strand
    findInSeq(seq);
    
    // Reverse complement strand (now correctly handles IUPAC codes)
    const revCompSeq = seq.reverseComplement();
    findInSeq(revCompSeq);
    
    return proteins;
  }
}
