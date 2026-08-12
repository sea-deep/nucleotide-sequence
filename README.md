# nucleotide-sequence 🧬

Hello! **nucleotide-sequence** is a little javascript library for Node.js and the Browser. It's meant to make sequence manipulation, translation, and alignment a bit easier and faster by using `Uint8Array`s under the hood. 

I originally wrote something like this years ago, and decided to modernize it for 2026. Hopefully you find it useful for your biology pipelines!

## 📦 Installation

```bash
npm install nucleotide-sequence
```

## 📖 Basic Usage

You can create sequences, manipulate them, and calculate basic properties.

```typescript
import { Seq } from 'nucleotide-sequence';

// Load a sequence
const mySeq = new Seq('DNA').read('ATGCCTGGATGC');

// Basic manipulation
const comp = mySeq.complement();
const rComp = mySeq.reverseComplement();
const gc = mySeq.gcContent(); // returns 0.5 (50%)

console.log(rComp.sequence()); // Prints the sequence string
```

## 🧬 Translation

You can easily translate sequences into amino acids using standard codon tables.

```typescript
import { Translation } from 'nucleotide-sequence';

const protein = new Translation(mySeq).translate();
console.log(protein.sequence());
```

## 🔍 Alignment (MatchMap)

If you need to find where a small sequence aligns within a larger sequence, you can use `MatchMap`. It handles wildcard matches (like `N`) as well.

```typescript
import { MatchMap } from 'nucleotide-sequence';

const query = new Seq().read('CCTG');
const map = new MatchMap(query, mySeq);
const result = map.initialize().best();

console.log(result.position); // Index where it best aligns
console.log(result.matches); // Number of matched nucleotides
```

## 🧵 Multithreading (ZeroWorker)

Genomic sequences can be quite large. If you are running an alignment in the browser and don't want to freeze the UI, you can offload it to a background thread using `ZeroWorker`. It tries to use `SharedArrayBuffer` or Transferable objects so it doesn't duplicate memory!

```typescript
import { ZeroWorker } from 'nucleotide-sequence';

const worker = new ZeroWorker();
const asyncResult = await worker.align(query, massiveGenome);
```

## 📚 Quick API Reference

### `Seq` Class
- **`new Seq(type)`**: Create a DNA or RNA sequence.
- **`.read(string)`**, **`.readFASTA(string)`**, **`.readFASTQ(string)`**: Load data.
- **`.sequence()`**: Get the string representation.
- **`.complement()`**, **`.reverseComplement()`**: Generate complements.
- **`.gcContent()`**: Calculate GC percentage.
- **`.splice(start, length)`**: Extract sub-regions instantly without copying memory.
- **`.kmers(k)`**: Extract all overlapping K-mers.
- **`.findCRISPRTargets(pam = 'NGG')`**: Scan for Cas9 PAM sites.

### `Translation` Engine
- **`new Translation(seq).translate()`**: Convert DNA/RNA to Amino Acids.

### `MatchMap`
- **`new MatchMap(query, reference).initialize().best()`**: Find the absolute best alignment position of a query inside a reference genome.

## 🐛 Feedback

This is a passion project, so there might be some edge cases I missed. If you run into any bugs or have suggestions, please open an issue on GitHub. I'd love to hear how you're using it!

## 📄 License
MIT License
