export type SeqType = 'DNA' | 'RNA';

/**
 * Nucleotide Sequence class optimized with Uint8Array for memory efficiency.
 */
export class Seq {
  public type: SeqType;
  private data: Uint8Array = new Uint8Array(0);

  constructor(type: SeqType = 'DNA') {
    this.type = type;
  }

  /**
   * Reads sequence data from a string.
   */
  public read(sequenceData: string): this {
    // Basic validation to clean up spaces and newlines if needed
    const cleanSeq = sequenceData.replace(/\s+/g, '').toUpperCase();
    const encoder = new TextEncoder();
    this.data = encoder.encode(cleanSeq);
    return this;
  }

  /**
   * Returns the size of the sequence in nucleotides.
   */
  public size(): number {
    return this.data.length;
  }

  /**
   * Returns the nucleotide sequence as a string.
   */
  public sequence(): string {
    const decoder = new TextDecoder();
    return decoder.decode(this.data);
  }

  /**
   * Returns a new Seq containing the complementary sequence.
   * A <-> T (or U in RNA)
   * C <-> G
   */
  public complement(): Seq {
    const compSeq = new Seq(this.type);
    const compData = new Uint8Array(this.data.length);

    for (let i = 0; i < this.data.length; i++) {
      const charCode = this.data[i];
      let compCode = charCode;

      switch (charCode) {
        case 65: // A
          compCode = this.type === 'DNA' ? 84 : 85; // T or U
          break;
        case 84: // T
        case 85: // U
          compCode = 65; // A
          break;
        case 67: // C
          compCode = 71; // G
          break;
        case 71: // G
          compCode = 67; // C
          break;
        case 77: // M (A/C) -> K (T/G)
          compCode = 75; 
          break;
        case 75: // K (G/T) -> M (C/A)
          compCode = 77; 
          break;
        case 82: // R (A/G) -> Y (T/C)
          compCode = 89; 
          break;
        case 89: // Y (C/T) -> R (G/A)
          compCode = 82; 
          break;
        case 66: // B (C/G/T) -> V (G/C/A)
          compCode = 86; 
          break;
        case 86: // V (A/C/G) -> B (T/G/C)
          compCode = 66; 
          break;
        case 68: // D (A/G/T) -> H (T/C/A)
          compCode = 72; 
          break;
        case 72: // H (A/C/T) -> D (T/G/A)
          compCode = 68; 
          break;
        // W, S, N are their own complements
        default:
          compCode = charCode;
      }
      compData[i] = compCode;
    }

    compSeq['data'] = compData;
    return compSeq;
  }

  /**
   * Transcribes DNA to RNA (T -> U).
   * Throws an error if called on an RNA sequence.
   */
  public transcribe(): Seq {
    if (this.type !== 'DNA') {
      throw new Error('Can only transcribe DNA to RNA');
    }
    const rna = new Seq('RNA');
    const rnaData = new Uint8Array(this.data.length);
    for (let i = 0; i < this.data.length; i++) {
      rnaData[i] = this.data[i] === 84 /* T */ ? 85 /* U */ : this.data[i];
    }
    rna['data'] = rnaData;
    return rna;
  }

  /**
   * Reverse Transcribes RNA to DNA (U -> T).
   */
  public reverseTranscribe(): Seq {
    if (this.type !== 'RNA') {
      throw new Error('Can only reverse transcribe RNA to DNA');
    }
    const dna = new Seq('DNA');
    const dnaData = new Uint8Array(this.data.length);
    for (let i = 0; i < this.data.length; i++) {
      dnaData[i] = this.data[i] === 85 /* U */ ? 84 /* T */ : this.data[i];
    }
    dna['data'] = dnaData;
    return dna;
  }

  /**
   * Calculates the GC content and other nucleotide fractions.
   */
  public fractionalContent(): Record<string, number> {
    const counts: Record<string, number> = { A: 0, C: 0, G: 0, T: 0, U: 0 };
    for (let i = 0; i < this.data.length; i++) {
      const char = String.fromCharCode(this.data[i]);
      if (counts[char] !== undefined) counts[char]++;
    }

    const total = this.data.length || 1;
    return {
      A: counts.A / total,
      C: counts.C / total,
      G: counts.G / total,
      T: counts.T / total,
      U: counts.U / total,
      GC: (counts.G + counts.C) / total,
    };
  }

  /**
   * NEW FEATURE: Melting Temperature (Tm) calculation using the Wallace Rule.
   * Useful for basic PCR primer design.
   */
  public meltingTemperature(): number {
    const counts = { A: 0, C: 0, G: 0, T: 0, U: 0 };
    for (let i = 0; i < this.data.length; i++) {
      const char = String.fromCharCode(this.data[i]);
      if (counts[char as keyof typeof counts] !== undefined) {
        counts[char as keyof typeof counts]++;
      }
    }
    
    // Wallace rule: Tm = 2°C(A+T) + 4°C(G+C)
    // Works well for primers 14-20 nucleotides long.
    const aTCount = counts.A + counts.T + counts.U;
    const gCCount = counts.G + counts.C;
    
    return (2 * aTCount) + (4 * gCCount);
  }

  /**
   * FASTA parser: Strips header lines starting with ">" and reads the raw sequence.
   */
  public readFASTA(fastaString: string): this {
    const lines = fastaString.split('\n');
    let sequence = '';
    for (const line of lines) {
      const trimmed = line.trim();
      if (!trimmed || trimmed.startsWith('>')) continue;
      sequence += trimmed;
    }
    return this.read(sequence);
  }

  /**
   * Calculates the approximate molecular weight (mass) in Daltons (g/mol).
   */
  public molecularWeight(): number {
    let mass = 0;
    for (let i = 0; i < this.data.length; i++) {
      const code = this.data[i];
      if (code === 65) mass += 313.2; // A
      else if (code === 84) mass += 304.2; // T
      else if (code === 67) mass += 289.2; // C
      else if (code === 71) mass += 329.2; // G
      else if (code === 85) mass += 306.2; // U
    }
    return mass;
  }

  /**
   * Removes or replaces existing nucleotides and/or adds new nucleotides in place.
   * Returns a NEW Seq object (immutable).
   */
  public splice(start: number, deleteCount: number = 0, insertSeq?: Seq): Seq {
    const sLen = this.data.length;
    let actualStart = start < 0 ? Math.max(sLen + start, 0) : Math.min(start, sLen);
    const actualDeleteCount = Math.min(Math.max(deleteCount, 0), sLen - actualStart);
    const insertLen = insertSeq ? insertSeq.size() : 0;
    
    const newData = new Uint8Array(sLen - actualDeleteCount + insertLen);
    
    // Copy left side
    newData.set(this.data.subarray(0, actualStart), 0);
    // Copy inserted sequence
    if (insertSeq) {
      newData.set(insertSeq['data'], actualStart);
    }
    // Copy right side
    newData.set(this.data.subarray(actualStart + actualDeleteCount), actualStart + insertLen);
    
    const newSeq = new Seq(this.type);
    newSeq['data'] = newData;
    return newSeq;
  }

  /**
   * FASTQ parser: Extracts only the raw sequence lines from a FASTQ format string.
   * Note: This discards the quality scores to save memory.
   */
  public readFASTQ(fastqString: string): this {
    const lines = fastqString.split('\n');
    let sequence = '';
    // In FASTQ, the 2nd line of every 4-line block is the sequence
    for (let i = 1; i < lines.length; i += 4) {
      sequence += lines[i].trim();
    }
    return this.read(sequence);
  }

  /**
   * Generates all overlapping K-mers of a given length `k`.
   * Extremely useful for modern alignment, Bloom filters, and ML models.
   */
  public kmers(k: number): string[] {
    const result: string[] = [];
    const seqStr = this.sequence();
    if (k > seqStr.length || k <= 0) return result;
    
    for (let i = 0; i <= seqStr.length - k; i++) {
      result.push(seqStr.slice(i, i + k));
    }
    return result;
  }

  /**
   * Scans the sequence for CRISPR PAM (Protospacer Adjacent Motif) target sites.
   * Defaults to 'NGG' for SpCas9. Returns a list of 0-indexed positions where the PAM starts.
   * Supports standard IUPAC degenerate codes (N, R, Y, W, S).
   */
  public findCRISPRTargets(pam: string = 'NGG'): number[] {
    // Convert IUPAC codes to Regex
    const regexStr = pam.toUpperCase()
      .replace(/N/g, '.')
      .replace(/R/g, '[AG]')
      .replace(/Y/g, '[CT]')
      .replace(/W/g, '[AT]')
      .replace(/S/g, '[GC]');
      
    // Use lookahead to find overlapping targets
    const regex = new RegExp(`(?=(${regexStr}))`, 'gi');
    const seqStr = this.sequence();
    const targets: number[] = [];
    
    let match;
    while ((match = regex.exec(seqStr)) !== null) {
      targets.push(match.index);
      // Advance manually to prevent infinite loop due to lookahead
      regex.lastIndex++;
    }
    return targets;
  }

  /**
   * Calculates the GC content (percentage of Guanine and Cytosine bases) in the sequence.
   * Returns a value between 0.0 and 1.0.
   */
  public gcContent(): number {
    let gcCount = 0;
    const len = this.data.length;
    for (let i = 0; i < len; i++) {
      const b = this.data[i];
      if (b === 71 /* G */ || b === 67 /* C */) {
        gcCount++;
      }
    }
    return len === 0 ? 0 : gcCount / len;
  }

  /**
   * Generates the reverse complement of the sequence.
   * Returns a new Seq object.
   */
  public reverseComplement(): Seq {
    const len = this.data.length;
    const revCompData = new Uint8Array(len);
    
    for (let i = 0; i < len; i++) {
      const b = this.data[len - 1 - i];
      let comp = b;
      
      // A (65) <-> T (84) / U (85)
      // C (67) <-> G (71)
      if (b === 65) comp = this.seqType === 'DNA' ? 84 : 85;
      else if (b === 84 || b === 85) comp = 65;
      else if (b === 67) comp = 71;
      else if (b === 71) comp = 67;
      
      revCompData[i] = comp;
    }

    const revSeq = new Seq(this.seqType);
    revSeq['data'] = revCompData;
    return revSeq;
  }
}
