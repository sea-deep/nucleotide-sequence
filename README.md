# nucleotide-sequence 🧬

[![NPM Version](https://img.shields.io/npm/v/nucleotide-sequence.svg)](https://www.npmjs.com/package/nucleotide-sequence)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Hello! **nucleotide-sequence** is a little javascript library for Node.js and the Browser. It's meant to make sequence manipulation, translation, and alignment a bit easier and faster by using `Uint8Array`s under the hood. 

I originally wrote something like this years ago, and decided to modernize it for 2026. Hopefully you find it useful for your biology pipelines!

---

## 📑 Table of Contents
- [Installation](#-installation)
- [Why Use This? (The Memory Problem)](#-why-use-this-the-memory-problem)
- [API Reference](#-api-reference)
  - [Sequence Engine (`Seq`)](#sequence-engine-seq)
  - [Translation Engine](#translation-engine)
  - [Alignment Mapping (`MatchMap`)](#alignment-mapping-matchmap)
- [License](#-license)

---

## 📦 Installation

```bash
npm install nucleotide-sequence
```

For browser environments (ESM):
```javascript
import { Seq, MatchMap, Translation, ZeroWorker } from 'nucleotide-sequence';
```

---

## 🤔 Why Use This? (The Memory Problem)

When dealing with gigabyte-scale genomic data in standard JavaScript, developers run into a catastrophic bottleneck: **Memory Bloat**.

Standard strings in JavaScript are UTF-16, meaning a 3GB genome takes up 6GB of RAM. If you try to pass that string to a Web Worker for background processing, the browser executes the **Structured Clone Algorithm**, copying it again. You instantly hit 12GB of RAM, and the Chrome tab crashes with an **"Aw, Snap! Out of Memory"** error.

**nucleotide-sequence** solves this by storing all biological data exclusively in flat, binary `Uint8Array` memory structures. 

---

## 📚 API Reference

### Sequence Engine (`Seq`)

The `Seq` class is the foundation of the library.

```typescript
import { Seq } from 'nucleotide-sequence';

const mySeq = new Seq('DNA').read('ATGCCTGGATGC');

const comp = mySeq.complement();
const rComp = mySeq.reverseComplement();
const gc = mySeq.gcContent(); // 0.5 (50%)
const rna = mySeq.transcribe();
```

**Advanced Methods:**
- `.readFASTQ(fastqData)`: Parses raw sequencer output while intelligently discarding bulky Phred quality scores.
- `.splice(start, length)`: Extracts sub-regions instantly using underlying ArrayBuffer slices (O(1) memory).
- `.kmers(k)`: Extracts overlapping K-mers.
- `.findCRISPRTargets(pam = 'NGG')`: Scans the sequence for CRISPR Cas9 targets.

### Translation Engine

Translate sequences into amino acids using standard biological codon tables.

```typescript
import { Translation } from 'nucleotide-sequence';

const protein = new Translation(mySeq).translate();
console.log(protein.sequence()); // Prints Amino Acid sequence
```

### Alignment Mapping (`MatchMap`)

Provides an exhaustive, sliding-window bitwise alignment engine to map small query sequences to massive reference genomes, fully supporting degenerate IUPAC nucleotides (like `N`).

```typescript
import { MatchMap } from 'nucleotide-sequence';

const query = new Seq().read('CCTG');
const map = new MatchMap(query, mySeq);
const result = map.initialize().best();

console.log(result.position); // Index where it best aligns
console.log(result.matches); // Number of matched nucleotides
```

---

## 🐛 Feedback

This is a passion project, so there might be some edge cases I missed. If you run into any bugs or have suggestions, please open an issue on GitHub. I'd love to hear how you're using it!

## 📄 License
MIT License
