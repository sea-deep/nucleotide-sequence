# nucleotide-sequence

[![NPM Version](https://img.shields.io/npm/v/nucleotide-sequence.svg)](https://www.npmjs.com/package/nucleotide-sequence)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**nucleotide-sequence** is a JavaScript/TypeScript library for Node.js and the Browser that provides functions for manipulating and analyzing DNA and RNA sequences. It uses `Uint8Array` to represent sequences internally.

## Installation

```bash
npm install nucleotide-sequence
```

## Basic Usage

```typescript
import { Seq, Translation, Alignment } from 'nucleotide-sequence';

// Initialize and read a sequence
const dna = new Seq('DNA').read('ATGCGTACGTTAG');

// Reverse Complement
const revComp = dna.reverseComplement();
console.log(revComp.sequence());

// Translation to Amino Acids (NCBI Table 1)
const protein = Translation.translate(dna);

// Pairwise Sequence Alignment (Smith-Waterman)
const ref = new Seq().read('ATGCGTACGT');
const result = Alignment.smithWaterman(dna, ref, { match: 2, mismatch: -1 });
console.log(`Alignment Score: ${result.score}`);
```

## API Reference

### `Seq`
The core class for wrapping and manipulating nucleotide sequences.
- **`read(sequence: string)`**: Parses a string into a `Uint8Array` sequence, ignoring whitespace.
- **`static readFASTA(content: string)`**: Parses a FASTA string and returns an array of `Seq` objects. Loads all records into memory; not suitable for gigabyte-scale genomic assemblies.
- **`static readFASTQ(content: string)`**: Parses a FASTQ string and returns an array of `Seq` objects.
- **`reverseComplement()`**: Returns a new `Seq` object containing the reverse complement, supporting IUPAC degenerate bases.
- **`kmers(k: number)`**: Returns a Generator yielding `Uint8Array` subarrays of length `k`. Yields strand-specific, overlapping substrings (not canonical k-mers).
- **`gcContent()`**: Computes the global GC percentage (0.0 to 1.0) based on exact G/C/S bases. Degenerate bases like 'N' are ignored.
- **`gcSkew(windowSize?: number)`**: Calculates `(G-C)/(G+C)` across sliding windows.
- **`hammingDistance(other: Seq)`**: Computes the Hamming distance between two sequences of equal length.
- **`meltingTemperatureNN(naConc?, primerConc?)`**: Computes primer Melting Temperature (Tm) using the SantaLucia (1998) Nearest-Neighbor parameters. Does not account for Mg²⁺ or sequence symmetry.
- **`molecularWeight({ phosphorylated?: boolean })`**: Computes mass using anhydrous IDT constants, applying a -61.96 Da correction to accurately model a standard linear synthetic oligo (5'-OH, 3'-OH).

### `Translation`
- **`static translate(seq: Seq, tableId?: number)`**: Translates a DNA/RNA sequence into an amino acid string using NCBI Translation Tables. Supports Standard (1), Vertebrate Mitochondrial (2), and Bacterial/Archaeal/Plant Plastid (11).
- **`static findOpenReadingFrames(seq: Seq, minLength?: number)`**: Scans all 6 reading frames to extract structural Open Reading Frames (ORFs) from start to stop codon. Does not predict biological translation viability (e.g., ignores Kozak/RBS sequences and eukaryotic introns).

### `Alignment`
- **`static smithWaterman(query: Seq, reference: Seq, options?: AlignmentOptions)`**: Performs local pairwise sequence alignment via dynamic programming.
- **`static needlemanWunsch(query: Seq, reference: Seq, options?: AlignmentOptions)`**: Performs global pairwise sequence alignment via dynamic programming.
*Note: Alignments utilize configurable integer match/mismatch scores and affine gap penalties. They do not implement evolutionary transition/transversion matrices like EDNAFULL.*

### `CrisprScoring`
- **`static findSpacers(seq: Seq, pam?: string, spacerLength?: number)`**: Identifies structural Protospacer Adjacent Motifs (PAMs) on both strands via regex. Does not evaluate chromatin accessibility.
- **`static calculateOnTargetScore(spacer: string)`**: Calculates an on-target efficiency score proxy based on a simplified positional weight matrix.
- **`static calculateCFDScore(guide: string, offTarget: string)`**: Calculates an off-target cutting frequency determination (CFD) score using static positional penalties.

### `SubstringSearch`
- **`constructor(query: Seq, reference: Seq)`**: Initializes an exact substring search tool.
- **`top(limit?: number)`**: Returns ungapped matches tolerant to `N` wildcards using an $O(M \times N)$ sliding window. Highly optimal for short amplicons and plasmids, but intractable for mapping short reads to whole genomes.

### `FormatSAM`
- **`static parse(samContent: string)`**: Parses Sequence Alignment/Map (SAM) text into an array of structured `SAMRecord` objects. Extracts standard fields but does not process bitwise FLAG semantics or CIGAR clipping operations.

### `Parallel`
*Requires the optional peer dependency `zeroworker`.*
- **`static align(query: Seq, references: Seq[], options?: AlignmentOptions)`**: Distributes pairwise alignments across a multithreaded Web Worker pool using `zeroworker`.
- **`static kmerCount(seq: Seq, k: number, chunks?: number)`**: Computes k-mer frequencies using a Web Worker pool.

## License
MIT License
