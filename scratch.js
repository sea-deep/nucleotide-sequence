const str = 'GGATCAGGAGGACTAGATGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTAACTGCAG';
const minLength = 15;
const startCodons = ['ATG', 'GTG', 'TTG'];
const stopCodons = ['TAA', 'TAG', 'TGA'];
const len = str.length;

console.log("length:", len);
console.log("ATG index:", str.indexOf('ATG'));
console.log("TAA index:", str.indexOf('TAA'));

for (let frame = 0; frame < 3; frame++) {
  let i = frame;
  while (i < len - 2) {
    const codon = str.substring(i, i + 3);
    if (startCodons.includes(codon)) {
      console.log(`Found start codon ${codon} at ${i}`);
      let j = i + 3;
      let foundStop = false;
      while (j < len - 2) {
        const nextCodon = str.substring(j, j + 3);
        if (stopCodons.includes(nextCodon)) {
          console.log(`Found stop codon ${nextCodon} at ${j}`);
          if (j - i >= minLength) {
            console.log("MATCH!");
          }
          i = j;
          foundStop = true;
          break;
        }
        j += 3;
      }
      if (!foundStop) {
        console.log("No stop found, breaking");
        break;
      }
    }
    i += 3;
  }
}
