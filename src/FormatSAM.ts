export interface SAMRecord {
  qname: string;
  flag: number;
  rname: string;
  pos: number;
  mapq: number;
  cigar: string;
  rnext: string;
  pnext: number;
  tlen: number;
  seq: string;
  qual: string;
  tags: Record<string, string | number>;
}

export class FormatSAM {
  /**
   * Parses a SAM file string into an array of SAMRecords.
   * Ignores header lines (starting with @).
   */
  public static parse(samContent: string): SAMRecord[] {
    const lines = samContent.split(/\r?\n/);
    const records: SAMRecord[] = [];

    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].trim();
      if (!line || line.startsWith('@')) continue;

      const fields = line.split('\t');
      if (fields.length < 11) {
        throw new Error(`Invalid SAM format at line ${i + 1}: expected at least 11 fields, found ${fields.length}`);
      }

      const record: SAMRecord = {
        qname: fields[0],
        flag: parseInt(fields[1], 10),
        rname: fields[2],
        pos: parseInt(fields[3], 10),
        mapq: parseInt(fields[4], 10),
        cigar: fields[5],
        rnext: fields[6],
        pnext: parseInt(fields[7], 10),
        tlen: parseInt(fields[8], 10),
        seq: fields[9],
        qual: fields[10],
        tags: {}
      };

      // Parse optional tags (TAG:TYPE:VALUE)
      for (let j = 11; j < fields.length; j++) {
        const tagParts = fields[j].split(':');
        if (tagParts.length >= 3) {
          const tag = tagParts[0];
          const type = tagParts[1];
          const value = tagParts.slice(2).join(':');

          if (type === 'i') {
            record.tags[tag] = parseInt(value, 10);
          } else if (type === 'f') {
            record.tags[tag] = parseFloat(value);
          } else {
            record.tags[tag] = value;
          }
        }
      }

      records.push(record);
    }

    return records;
  }
}
