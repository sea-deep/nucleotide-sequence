import { describe, it, expect } from 'vitest';
import { Seq } from '../src/Seq';
import { GenePrediction } from '../src/GenePrediction';

describe('GenePrediction', () => {
  it('should predict genes based on ORFs and codon bias', () => {
    const geneBody = 'CCG'.repeat(20); // 20 CCG codons
    const seqStr = 'GGATCAGGAGGACTAGATG' + geneBody + 'TAACTGCAG';
    
    const seq = new Seq().read(seqStr);
    const predictor = new GenePrediction(seq);
    
    const models = predictor.predict(5);
    
    expect(models.length).toBeGreaterThan(0);
    expect(models[0].strand).toBe('+');
    expect(models[0].score).toBeGreaterThan(0);
    // Should have detected the RBS
    expect(seqStr.substring(models[0].start, models[0].end).startsWith('ATG')).toBe(true);
  });
});
