export type SeqType = 'DNA' | 'RNA';

// Valid IUPAC nucleotide codes (uppercase ASCII)
const VALID_IUPAC = new Uint8Array(256);
for (const c of 'ACGTUMRWSYKVHDBN') {
  VALID_IUPAC[c.charCodeAt(0)] = 1;
}

// Complement lookup table (ASCII to ASCII)
const COMPLEMENT_TABLE = new Uint8Array(256);
// Initialize with self (fallback)
for (let i = 0; i < 256; i++) COMPLEMENT_TABLE[i] = i;
// Standard
COMPLEMENT_TABLE[65] = 84; // A -> T
COMPLEMENT_TABLE[84] = 65; // T -> A
COMPLEMENT_TABLE[85] = 65; // U -> A
COMPLEMENT_TABLE[67] = 71; // C -> G
COMPLEMENT_TABLE[71] = 67; // G -> C
// Degenerate
COMPLEMENT_TABLE[77] = 75; // M -> K
COMPLEMENT_TABLE[75] = 77; // K -> M
COMPLEMENT_TABLE[82] = 89; // R -> Y
COMPLEMENT_TABLE[89] = 82; // Y -> R
COMPLEMENT_TABLE[87] = 87; // W -> W
COMPLEMENT_TABLE[83] = 83; // S -> S
COMPLEMENT_TABLE[86] = 66; // V -> B
COMPLEMENT_TABLE[66] = 86; // B -> V
COMPLEMENT_TABLE[72] = 68; // H -> D
COMPLEMENT_TABLE[68] = 72; // D -> H
COMPLEMENT_TABLE[78] = 78; // N -> N

// IS_GC lookup table
const IS_GC = new Uint8Array(256);
IS_GC[67] = 1; // C
IS_GC[71] = 1; // G
IS_GC[83] = 1; // S (C/G)

/**
 * Nucleotide Sequence class using Uint8Array for storage.
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
    const cleanSeq = sequenceData.replace(/\s+/g, '').toUpperCase();
    const encoder = new TextEncoder();
    const encoded = encoder.encode(cleanSeq);

    // Validate every character against the IUPAC alphabet
    for (let i = 0; i < encoded.length; i++) {
      if (!VALID_IUPAC[encoded[i]]) {
        const char = String.fromCharCode(encoded[i]);
        throw new TypeError(
          `Invalid character '${char}' at position ${i}. Expected IUPAC nucleotide code (ACGTUMRWSYKVHDBN).`
        );
      }
    }

    this.data = encoded;
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
   * Generates the complementary sequence.
   * Uses a high-performance lookup table.
   * Returns a new Seq object.
   */
  public complement(): Seq {
    const len = this.data.length;
    const compData = new Uint8Array(len);
    
    for (let i = 0; i < len; i++) {
      let comp = COMPLEMENT_TABLE[this.data[i]];
      // Handle RNA specific T/U swap if necessary
      if (this.type === 'DNA' && comp === 85) comp = 84; // U -> T
      else if (this.type === 'RNA' && comp === 84) comp = 85; // T -> U
      compData[i] = comp;
    }

    const compSeq = new Seq(this.type);
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
    let a = 0, c = 0, g = 0, t = 0, u = 0;
    for (let i = 0; i < this.data.length; i++) {
      const b = this.data[i];
      if (b === 65) a++;      // A
      else if (b === 67) c++; // C
      else if (b === 71) g++; // G
      else if (b === 84) t++; // T
      else if (b === 85) u++; // U
    }

    const total = this.data.length || 1;
    return {
      A: a / total,
      C: c / total,
      G: g / total,
      T: t / total,
      U: u / total,
      GC: (g + c) / total,
    };
  }

  /**
   * Melting Temperature (Tm) estimation using the Wallace Rule.
   * Valid only for oligonucleotides of 14-20 bases (Wallace et al. 1979).
   * For longer sequences or precise Tm, use nearest-neighbor thermodynamics.
   *
   * @throws {RangeError} If the sequence length is outside the valid range.
   */
  public meltingTemperature(): number {
    const len = this.data.length;
    if (len < 14 || len > 20) {
      throw new RangeError(
        `Wallace Rule is valid for 14-20nt oligomers. This sequence is ${len}nt. ` +
        `For longer sequences, use nearest-neighbor thermodynamics.`
      );
    }

    let at = 0, gc = 0;
    for (let i = 0; i < len; i++) {
      const b = this.data[i];
      if (b === 65 || b === 84 || b === 85) at++; // A, T, U
      else if (b === 71 || b === 67) gc++;        // G, C
    }
    
    return (2 * at) + (4 * gc);
  }


  /**
   * Calculates the approximate molecular weight (mass) in Daltons (g/mol).
   * 
   * By default, this models a standard synthetic oligonucleotide (5'-OH, 3'-OH).
   * To model a 5'-phosphorylated oligo, pass `{ phosphorylated: true }`.
   * 
   * Note: This calculates the single-stranded mass using standard IDT constants.
   */
  public molecularWeight(options: { phosphorylated?: boolean } = {}): number {
    let mass = 0;
    for (let i = 0; i < this.data.length; i++) {
      const code = this.data[i];
      if (code === 65) mass += 313.21; // A (dAMP residue)
      else if (code === 84) mass += 304.20; // T (dTMP residue)
      else if (code === 67) mass += 289.18; // C (dCMP residue)
      else if (code === 71) mass += 329.21; // G (dGMP residue)
      else if (code === 85) mass += 306.17; // U (UMP residue)
    }

    if (this.data.length === 0) return 0;

    // The residue masses above inherently contain one phosphate per base.
    // A standard linear oligo (5'-OH, 3'-OH) has one LESS phosphate than bases,
    // plus two terminal hydroxyl/hydrogen atoms.
    // Correction = + H2O (18.02) - PO3 (79.98) = -61.96
    if (options.phosphorylated) {
      // 5'-PO4, 3'-OH: Just add the H2O to cap the ends.
      return mass + 18.02;
    } else {
      // 5'-OH, 3'-OH: Standard IDT formula
      return mass - 61.96;
    }
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
   * Parses a multi-record FASTA string.
   * Returns an array of objects containing the ID, description, and the parsed Seq object.
   * 
   * WARNING: This loads all sequences into an array of objects in memory. 
   * The memory footprint scales linearly with the sequence size. This is optimal 
   * for plasmids, amplicons, and viral databases, but unsuitable for 
   * gigabyte-scale genomic assemblies.
   */
  public static readFASTA(data: string, type: SeqType = 'DNA'): { id: string, description: string, seq: Seq }[] {
    const records: { id: string, description: string, seq: Seq }[] = [];
    const blocks = data.split(/(?:^|\n)>/);
    for (const block of blocks) {
      if (!block.trim()) continue;
      
      const lines = block.split('\n');
      const header = lines[0].trim();
      const id = header.split(/\s+/)[0];
      const description = header.substring(id.length).trim();
      
      const sequence = lines.slice(1).join('');
      const seq = new Seq(type).read(sequence);
      
      records.push({ id, description, seq });
    }
    return records;
  }

  /**
   * Parses a FASTQ string with quality scores.
   * Returns an array of FastqRecord objects.
   */
  public static readFASTQ(data: string, type: SeqType = 'DNA'): { id: string, seq: Seq, quality: Uint8Array }[] {
    const records: { id: string, seq: Seq, quality: Uint8Array }[] = [];
    const lines = data.split('\n');
    
    for (let i = 0; i < lines.length; i += 4) {
      if (!lines[i] || lines[i][0] !== '@') continue;
      
      const id = lines[i].substring(1).split(/\s+/)[0];
      const seqStr = lines[i + 1].trim();
      // Skip lines[i+2] which is just '+'
      const qualStr = lines[i + 3] ? lines[i + 3].trim() : '';
      
      if (seqStr.length !== qualStr.length) {
        throw new Error(`FASTQ format error: Sequence and quality score lengths differ for record ${id}`);
      }
      
      const seq = new Seq(type).read(seqStr);
      const encoder = new TextEncoder();
      const quality = encoder.encode(qualStr);
      
      records.push({ id, seq, quality });
    }
    return records;
  }

  /**
   * Computes the GC skew in sliding windows.
   * GC skew is (G-C)/(G+C) and is used to find origin/terminus of replication.
   */
  public gcSkew(windowSize: number): { position: number, skew: number }[] {
    if (windowSize <= 0 || windowSize > this.data.length) {
      throw new RangeError('Invalid window size for GC skew calculation');
    }
    
    const results: { position: number, skew: number }[] = [];
    for (let i = 0; i <= this.data.length - windowSize; i += windowSize) {
      let g = 0;
      let c = 0;
      for (let j = 0; j < windowSize; j++) {
        const b = this.data[i + j];
        if (b === 71) g++;
        else if (b === 67) c++;
      }
      const sum = g + c;
      const skew = sum === 0 ? 0 : (g - c) / sum;
      results.push({ position: i, skew });
    }
    return results;
  }

  /**
   * Calculates the Hamming distance between this sequence and another.
   * @throws {Error} if sequences are of unequal length.
   */
  public hammingDistance(other: Seq): number {
    if (this.data.length !== other.data.length) {
      throw new Error('Hamming distance requires sequences of equal length');
    }
    
    let distance = 0;
    const len = this.data.length;
    for (let i = 0; i < len; i++) {
      if (this.data[i] !== other.data[i]) {
        distance++;
      }
    }
    return distance;
  }

  /**
   * Calculates nearest-neighbor melting temperature.
   * Based on SantaLucia (1998) unified thermodynamic parameters.
   */
  public meltingTemperatureNN(naConc: number = 0.05, primerConc: number = 0.00000025): number {
    const len = this.data.length;
    if (len < 2) return 0;
    
    // Simplified unified NN parameters for DNA (SantaLucia 1998)
    // Values are roughly approximated for this implementation
    // {dH (kcal/mol), dS (cal/K*mol)}
    const nnParams: Record<string, {dH: number, dS: number}> = {
      'AA': {dH: -7.9, dS: -22.2}, 'TT': {dH: -7.9, dS: -22.2},
      'AT': {dH: -7.2, dS: -20.4},
      'TA': {dH: -7.2, dS: -21.3},
      'CA': {dH: -8.5, dS: -22.7}, 'TG': {dH: -8.5, dS: -22.7},
      'GT': {dH: -8.4, dS: -22.4}, 'AC': {dH: -8.4, dS: -22.4},
      'CT': {dH: -7.8, dS: -21.0}, 'AG': {dH: -7.8, dS: -21.0},
      'GA': {dH: -8.2, dS: -22.2}, 'TC': {dH: -8.2, dS: -22.2},
      'CG': {dH: -10.6, dS: -27.2},
      'GC': {dH: -9.8, dS: -24.4},
      'GG': {dH: -8.0, dS: -19.9}, 'CC': {dH: -8.0, dS: -19.9},
    };
    
    let sum_dH = 0; // cal/mol
    let sum_dS = 0; // cal/K*mol
    
    // Initiation parameters
    const init_dH = 0.2; // kcal/mol
    const init_dS = -5.7; // cal/K*mol
    sum_dH += init_dH * 1000;
    sum_dS += init_dS;
    
    // Convert array back to string just for nearest neighbor pairs
    const seqStr = this.sequence().toUpperCase();
    for (let i = 0; i < len - 1; i++) {
      const pair = seqStr.substring(i, i + 2);
      if (nnParams[pair]) {
        sum_dH += nnParams[pair].dH * 1000;
        sum_dS += nnParams[pair].dS;
      }
    }
    
    // Gas constant R = 1.987 cal/K*mol
    const R = 1.987;
    
    // Tm = (dH / (dS + R * ln(C))) - 273.15 + 16.6 * log10([Na+])
    // C is the total primer concentration
    const tm = (sum_dH / (sum_dS + R * Math.log(primerConc))) - 273.15 + 16.6 * Math.log10(naConc);
    return Math.max(0, tm); // Ensure no negative Tm
  }


  /**
   * Extracts all overlapping k-mers of a given size.
   * Returns a Generator of Uint8Array views (zero-copy subarrays).
   * @param k The size of the k-mer.
   */
  public *kmers(k: number): Generator<Uint8Array> {
    if (k <= 0 || k > this.data.length) {
      throw new Error('Invalid k-mer size');
    }
    
    for (let i = 0; i <= this.data.length - k; i++) {
      yield this.data.subarray(i, i + k);
    }
  }

  /**
   * Scans the sequence for structural PAM (Protospacer Adjacent Motif) sites.
   * Defaults to 'NGG' (SpCas9). Returns a list of 0-indexed positions where the PAM starts.
   * Supports standard IUPAC degenerate codes (N, R, Y, W, S).
   *
   * Note: This finds structural PAM motif occurrences only. It does not predict CRISPR
   * cleavage efficiency, chromatin accessibility, or off-target thermodynamics.
   */
  public findPAMSites(pam: string = 'NGG'): number[] {
    const regexStr = pam.toUpperCase()
      .replace(/N/g, '.')
      .replace(/R/g, '[AG]')
      .replace(/Y/g, '[CT]')
      .replace(/W/g, '[AT]')
      .replace(/S/g, '[GC]');
      
    const regex = new RegExp(`(?=(${regexStr}))`, 'gi');
    const seqStr = this.sequence();
    const sites: number[] = [];
    
    let match;
    while ((match = regex.exec(seqStr)) !== null) {
      sites.push(match.index);
      regex.lastIndex++;
    }
    return sites;
  }

  /**
   * Calculates the GC content of the sequence.
   * Uses lookup table for performance.
   * @returns A float between 0 and 1.
   */
  public gcContent(): number {
    let gcCount = 0;
    const len = this.data.length;
    for (let i = 0; i < len; i++) {
      if (IS_GC[this.data[i]]) {
        gcCount++;
      }
    }
    return len === 0 ? 0 : gcCount / len;
  }

  /**
   * Generates the reverse complement of the sequence.
   * Handles all IUPAC degenerate nucleotide codes.
   * Returns a new Seq object.
   */
  public reverseComplement(): Seq {
    const comp = this.complement();
    const len = comp['data'].length;
    const reversed = new Uint8Array(len);
    for (let i = 0; i < len; i++) {
      reversed[i] = comp['data'][len - 1 - i];
    }
    const revSeq = new Seq(this.type);
    revSeq['data'] = reversed;
    return revSeq;
  }
}
