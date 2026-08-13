# nucleotide-sequence 🧬

[![NPM Version](https://img.shields.io/npm/v/nucleotide-sequence.svg)](https://www.npmjs.com/package/nucleotide-sequence)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**nucleotide-sequence** is a modern, high-performance JavaScript/TypeScript library for working with DNA and RNA sequences. It provides standard bioinformatics operations including alignment, translation, CRISPR targeting, and thermodynamic properties.

Built for both Node.js and the browser, the library operates internally on zero-copy `Uint8Array` data structures, minimizing V8 heap allocation and garbage collection pauses to process genomic data at exceptional speeds.

**Note:** This is a utility library suitable for analyzing plasmids, viral genomes, amplicons, and individual genes. It is *not* designed for massive genome-scale indexing or millions-of-reads processing pipelines.

## 📦 Installation

```bash
npm install nucleotide-sequence
```

## 📖 Basic Usage

```typescript
import { Seq, Translation, Alignment } from 'nucleotide-sequence';

// 1. Initialize and read a sequence
const dna = new Seq('DNA').read('ATGCGTACGTTAG');

// 2. Fast Reverse Complement
const revComp = dna.reverseComplement();
console.log(revComp.sequence());

// 3. Translation to Amino Acids (NCBI Table 1)
const protein = Translation.translate(dna);

// 4. Pairwise Sequence Alignment (Smith-Waterman)
const ref = new Seq().read('ATGCGTACGT');
const result = Alignment.smithWaterman(dna, ref, { match: 2, mismatch: -1 });
console.log(`Alignment Score: ${result.score}`);
```

## 📚 API Reference

### `Seq`
The core class for wrapping and manipulating nucleotide sequences.
- **`read(sequence: string)`**: Parses a string into a `Uint8Array` sequence, stripping whitespace automatically.
- **`static readFASTA(content: string)`**: Parses a FASTA file and returns an array of `Seq` objects.
- **`static readFASTQ(content: string)`**: Parses a FASTQ file and returns an array of `Seq` objects.
- **`reverseComplement()`**: Returns a new `Seq` object containing the reverse complement, supporting all IUPAC degenerate bases.
- **`kmers(k: number)`**: Returns a Generator yielding zero-copy `Uint8Array` subarrays for highly efficient k-mer iteration.
- **`gcContent()`**: Computes the global GC percentage using a fast byte-lookup table.
- **`gcSkew(windowSize?: number)`**: Calculates `(G-C)/(G+C)` across sliding windows to identify origins of replication.
- **`hammingDistance(other: Seq)`**: Computes the Hamming distance between two sequences.
- **`meltingTemperatureNN(dnaConcentration?, saltConcentration?)`**: Computes standard primer Melting Temperature (Tm) using the SantaLucia (1998) Nearest-Neighbor thermodynamic parameters.

### `Translation`
- **`static translate(seq: Seq, tableId?: number)`**: Translates a DNA/RNA sequence into an amino acid string using NCBI Translation Tables. Supports Standard (1), Vertebrate Mitochondrial (2), and Bacterial/Archaeal/Plant Plastid (11) by default.
- **`static findOpenReadingFrames(seq: Seq, minLength?: number)`**: Scans all 6 reading frames to extract structural Open Reading Frames (ORFs).

### `Alignment`
- **`static smithWaterman(query: Seq, reference: Seq, options?: AlignmentOptions)`**: Performs local pairwise sequence alignment via dynamic programming and affine gap penalties.
- **`static needlemanWunsch(query: Seq, reference: Seq, options?: AlignmentOptions)`**: Performs global pairwise sequence alignment via dynamic programming and affine gap penalties.

### `CrisprScoring`
- **`static extractSpacers(seq: Seq, pam?: string, spacerLength?: number)`**: Identifies Protospacer Adjacent Motifs (PAMs) on both the forward and reverse strands, extracting the adjacent spacer sequences and calculating their coordinates and GC bounds.

### `SubstringSearch`
- **`constructor(query: Seq, reference: Seq)`**: Initializes an exact substring search tool.
- **`top(limit?: number)`**: Returns the top ungapped alignments, tolerant to `N` wildcards.

### `Parallel`
*Requires the optional peer dependency `zeroworker`.*
- **`static align(query: Seq, references: Seq[], options?: AlignmentOptions)`**: Distributes pairwise alignments across a multithreaded Web Worker pool.
- **`static kmerCount(seq: Seq, k: number)`**: Performs parallelized k-mer counting.

## 🐛 Feedback & Contributions
This project is open-source. If you find bugs or have feature requests, please open an issue on GitHub.

## 📄 License
MIT License
