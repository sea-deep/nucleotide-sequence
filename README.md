# nucleotide-sequence 🧬

**nucleotide-sequence** is a high-performance bioinformatics library written in TypeScript, providing sequence manipulation, analysis, and bitwise alignment mapping for Node.js and the browser.

Built to handle large-scale genomic datasets, it utilizes `Uint8Array` memory structures to minimize overhead compared to traditional string-based approaches. The library features built-in support for FASTA/FASTQ parsing, CRISPR PAM scanning, K-mer extraction, and exhaustive degenerate nucleotide alignment.

## 🚀 Features

- **Blazing Fast**: Uses Typed Arrays (`Uint8Array`) to minimize memory footprint and garbage collection overhead.
- **MatchMap Algorithm**: Exhaustively scans genomes using an optimized sliding window, running at millions of nucleotide comparisons per second.
- **Modern Standards**: Seamlessly extract **K-mers**, parse raw **FASTQ** sequencer outputs, and locate **CRISPR** editable targets (e.g. `NGG`).
- **Translation Engine**: Instantly translate DNA/RNA into Amino Acids, or scan all 6 reading frames for hidden **Open Reading Frames (ORFs)**.
- **Degenerate Nucleotides**: Full support for standard IUPAC codes (`N`, `R`, `Y`, `W`, `S`, etc.).
- **Zero Dependencies**: A lightweight, standalone library built for the modern web and server.

## 📦 Installation

**nucleotide-sequence** is available on NPM:

```bash
npm install nucleotide-sequence
```

For browser environments, you can import it directly via ESM:

```javascript
import { Seq, MatchMap, Translation } from 'nucleotide-sequence';
```

---

## 📖 Quick Usage

```typescript
import { Seq, MatchMap, Translation } from 'nucleotide-sequence';

// 1. Create and manipulate a sequence
const seq = new Seq().read('ATGCATGC');
const rna = seq.transcribe();
console.log(rna.sequence()); // "AUGCAUGC"

// 2. Parse a raw FASTQ file
const fastq = `@SEQ_ID\nGATTTGGGGTTCAA\n+\n!''*((((***+))`;
const rawSeq = new Seq().readFASTQ(fastq);

// 3. Find CRISPR Targets
const crisprTargets = rawSeq.findCRISPRTargets('NGG');
console.log(`Found SpCas9 cut sites at indices:`, crisprTargets);

// 4. Extract K-mers for ML/Bloom Filters
const kmers = seq.kmers(4); // ["ATGC", "TGCA", "GCAT", ...]

// 5. High-Performance Alignment Mapping
const reference = new Seq().read('ATGCCTGGATGC');
const query = new Seq().read('ATGC');
const map = new MatchMap(query, reference).initialize();

// Get the absolute best match
const bestMatch = map.best();
console.log(`Matched at index ${bestMatch?.position} with ${bestMatch?.matches} exact matches!`);
```

---

## 📚 Library Reference

### `Seq` Class

#### `new Seq(type?: 'DNA' | 'RNA')`
Constructs a new sequence. Defaults to `'DNA'`.

#### `.read(data: string): this`
Reads string data into the optimized `Uint8Array`. Automatically ignores whitespace and newlines.

#### `.readFASTA(fastaData: string): this`
Parses a standard `.fasta` file, stripping headers and extracting the raw sequence.

#### `.readFASTQ(fastqData: string): this`
Parses raw output from modern sequencers, extracting the sequence while discarding the quality scores for memory efficiency.

#### `.sequence(): string`
Returns the nucleotide sequence as a standard string.

#### `.complement(): Seq`
Returns a new `Seq` object containing the 5' -> 3' complement. Properly handles RNA and DNA degenerate nucleotides.

#### `.transcribe(): Seq`
Returns a new `Seq` of type `'RNA'` (converts `T` to `U`).

#### `.reverseTranscribe(): Seq`
Returns a new `Seq` of type `'DNA'` (converts `U` to `T`).

#### `.splice(start: number, deleteCount: number = 0, insertSeq?: Seq): Seq`
A highly efficient typed array mutation function. Removes/replaces nucleotides and/or adds new nucleotides in place, returning a new immutable `Seq`.

#### `.gcContent(): number`
Returns the percentage of G and C nucleotides in the sequence (0.0 to 1.0).

#### `.meltingTemperature(): number`
Calculates the Wallace Rule melting temperature (Tm) for PCR primers.

#### `.molecularWeight(): number`
Calculates the physical mass of the sequence in Daltons (g/mol).

#### `.kmers(k: number): string[]`
Extracts all overlapping K-mers of length `k`. Essential for ML tokenization and alignment building.

#### `.findCRISPRTargets(pam: string = 'NGG'): number[]`
Scans the genome for CRISPR PAM (Protospacer Adjacent Motif) target sites. Returns an array of 0-indexed positions where the PAM starts. Supports standard IUPAC degenerate codes.

---

### `Translation` Engine

#### `Translation.translate(seq: Seq): string`
Translates a DNA or RNA sequence into its corresponding Amino Acid string using the standard biological codon table. Unrecognized/degenerate codons are returned as `?`.

#### `Translation.translateFrame(seq: Seq, frame: number): string`
Translates the sequence starting from a specific reading frame (0, 1, or 2).

#### `Translation.findOpenReadingFrames(seq: Seq): string[]`
Scans all 6 reading frames (3 forward, 3 reverse complement) to extract valid hidden proteins. It searches for sequences that start with Methionine (`M`) and end with a Stop Codon (`*`).

---

### `MatchMap` Algorithm

The MatchMap algorithm provides a highly optimized sliding-window exhaustive bitwise alignment.

#### `new MatchMap(query: Seq, reference: Seq)`
Initializes the engine to search for the `query` inside the massive `reference` genome.

#### `.initialize(): this`
Executes the alignment map.

#### `.best(): MatchResult | null`
Returns the single best alignment match.

#### `.top(maxResults: number): MatchResult[]`
Returns the top `N` alignments, sorted by the highest number of matches.

---

### `MatchResult` Class
Returned by the `MatchMap` queries.

- `.position`: The 0-indexed alignment position in the reference sequence.
- `.matches`: The number of successful matches (including `N` wildcards).
- `.alignment()`: Returns a new `Seq` containing the exactly matched portion of the reference sequence.

---

## 📄 License
MIT License
