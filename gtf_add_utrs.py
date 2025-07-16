#!/usr/bin/env python3
"""
Script to add UTRs to GTF annotations where exon coordinates match CDS coordinates.
Adds configurable length 5' UTR and 3' UTR by extending terminal exons.
Default: 150bp 5' UTR and 250bp 3' UTR.
"""

import sys
import re
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple

class GTFRecord:
    def __init__(self, line: str):
        fields = line.strip().split('\t')
        self.seqname = fields[0]
        self.source = fields[1]
        self.feature = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        self.attribute = fields[8]
        
        # Parse attributes
        self.attributes = {}
        for attr in self.attribute.split(';'):
            attr = attr.strip()
            if attr:
                match = re.match(r'(\w+)\s+"([^"]+)"', attr)
                if match:
                    self.attributes[match.group(1)] = match.group(2)
    
    def get_transcript_id(self) -> str:
        return self.attributes.get('transcript_id', '')
    
    def get_gene_id(self) -> str:
        return self.attributes.get('gene_id', '')
    
    def to_line(self) -> str:
        return f"{self.seqname}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.frame}\t{self.attribute}"

def parse_gtf(filename: str) -> List[GTFRecord]:
    """Parse GTF file and return list of GTFRecord objects."""
    records = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                records.append(GTFRecord(line))
    return records

def parse_chromosome_sizes(filename: str) -> Dict[str, int]:
    """Parse chromosome sizes file and return dictionary of chromosome -> size."""
    chrom_sizes = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom_name = parts[0]
                    try:
                        chrom_size = int(parts[1])
                        chrom_sizes[chrom_name] = chrom_size
                    except ValueError:
                        print(f"Warning: Invalid chromosome size for {chrom_name}: {parts[1]}")
    return chrom_sizes

def group_by_transcript(records: List[GTFRecord]) -> Dict[str, List[GTFRecord]]:
    """Group GTF records by transcript ID."""
    transcript_groups = defaultdict(list)
    for record in records:
        transcript_id = record.get_transcript_id()
        if transcript_id:
            transcript_groups[transcript_id].append(record)
    return transcript_groups

def ranges_match(exon_ranges: List[Tuple[int, int]], cds_ranges: List[Tuple[int, int]]) -> bool:
    """Check if exon ranges exactly match CDS ranges."""
    if len(exon_ranges) != len(cds_ranges):
        return False
    
    # Sort ranges by start position
    exon_ranges_sorted = sorted(exon_ranges)
    cds_ranges_sorted = sorted(cds_ranges)
    
    return exon_ranges_sorted == cds_ranges_sorted

def create_utr_record(template: GTFRecord, feature_type: str, start: int, end: int) -> GTFRecord:
    """Create a new UTR record based on a template record."""
    utr_record = GTFRecord(template.to_line())
    utr_record.feature = feature_type
    utr_record.start = start
    utr_record.end = end
    utr_record.frame = '.'  # UTRs don't have frame
    
    # Rebuild attribute string
    attr_parts = []
    for key, value in utr_record.attributes.items():
        attr_parts.append(f'{key} "{value}"')
    utr_record.attribute = '; '.join(attr_parts) + ';'
    
    return utr_record

def process_transcript(transcript_records: List[GTFRecord], utr_5_length: int = 150, utr_3_length: int = 250, chrom_sizes: Dict[str, int] = None) -> List[GTFRecord]:
    """Process a single transcript and add UTRs if exon ranges match CDS ranges"""
    
    # Separate records by feature type
    transcript_record = None
    exon_records = []
    cds_records = []
    other_records = []
    
    for record in transcript_records:
        if record.feature == 'transcript':
            transcript_record = record
        elif record.feature == 'exon':
            exon_records.append(record)
        elif record.feature == 'CDS':
            cds_records.append(record)
        else:
            other_records.append(record)
    
    # Check if we have the required features
    if not transcript_record or not exon_records or not cds_records:
        return transcript_records
    
    # Get ranges
    exon_ranges = [(r.start, r.end) for r in exon_records]
    cds_ranges = [(r.start, r.end) for r in cds_records]
    
    # Check if exon ranges match CDS ranges
    if not ranges_match(exon_ranges, cds_ranges):
        return transcript_records
    
    # Sort exon records by start position
    exon_records_sorted = sorted(exon_records, key=lambda x: x.start)
    strand = transcript_record.strand
    chromosome = transcript_record.seqname
    
    # Get chromosome size if available
    chrom_size = chrom_sizes.get(chromosome) if chrom_sizes else None
    
    # Determine terminal exons based on strand
    if strand == '+':
        first_exon = exon_records_sorted[0]
        last_exon = exon_records_sorted[-1]
        
        # Calculate 5' UTR extension (upstream) - ensure it doesn't go below 1
        original_5_start = first_exon.start
        desired_5_start = original_5_start - utr_5_length
        actual_5_start = max(1, desired_5_start)  # GTF is 1-based
        actual_5_length = original_5_start - actual_5_start
        
        # Calculate 3' UTR extension (downstream) - ensure it doesn't exceed chromosome size
        original_3_end = last_exon.end
        desired_3_end = original_3_end + utr_3_length
        if chrom_size:
            actual_3_end = min(chrom_size, desired_3_end)
        else:
            actual_3_end = desired_3_end
        actual_3_length = actual_3_end - original_3_end
        
        # UTR coordinates
        utr_5_start = actual_5_start
        utr_5_end = original_5_start - 1
        utr_3_start = original_3_end + 1
        utr_3_end = actual_3_end
        
        # Update exon coordinates
        first_exon.start = actual_5_start
        last_exon.end = actual_3_end
        
        # Report if UTRs were truncated
        if actual_5_length < utr_5_length:
            print(f"    Warning: 5' UTR truncated to {actual_5_length}bp (requested {utr_5_length}bp) for {transcript_record.get_transcript_id()}")
        if chrom_size and actual_3_length < utr_3_length:
            print(f"    Warning: 3' UTR truncated to {actual_3_length}bp (requested {utr_3_length}bp) for {transcript_record.get_transcript_id()}")
        
    else:  # strand == '-'
        first_exon = exon_records_sorted[0]
        last_exon = exon_records_sorted[-1]
        
        # For - strand, 5' UTR is downstream (after last exon), 3' UTR is upstream (before first exon)
        original_3_start = first_exon.start
        original_5_end = last_exon.end
        
        # Calculate 3' UTR extension (upstream) - ensure it doesn't go below 1
        desired_3_start = original_3_start - utr_3_length
        actual_3_start = max(1, desired_3_start)
        actual_3_length = original_3_start - actual_3_start
        
        # Calculate 5' UTR extension (downstream) - ensure it doesn't exceed chromosome size
        desired_5_end = original_5_end + utr_5_length
        if chrom_size:
            actual_5_end = min(chrom_size, desired_5_end)
        else:
            actual_5_end = desired_5_end
        actual_5_length = actual_5_end - original_5_end
        
        # UTR coordinates
        utr_3_start = actual_3_start
        utr_3_end = original_3_start - 1
        utr_5_start = original_5_end + 1
        utr_5_end = actual_5_end
        
        # Update exon coordinates
        first_exon.start = actual_3_start
        last_exon.end = actual_5_end
        
        # Report if UTRs were truncated
        if actual_3_length < utr_3_length:
            print(f"    Warning: 3' UTR truncated to {actual_3_length}bp (requested {utr_3_length}bp) for {transcript_record.get_transcript_id()}")
        if chrom_size and actual_5_length < utr_5_length:
            print(f"    Warning: 5' UTR truncated to {actual_5_length}bp (requested {utr_5_length}bp) for {transcript_record.get_transcript_id()}")
    
    # update transcript coordinates
    all_starts = [r.start for r in exon_records_sorted]
    all_ends = [r.end for r in exon_records_sorted]
    transcript_record.start = min(all_starts)
    transcript_record.end = max(all_ends)
    
    # create utr records only if they have positive length
    new_records = []
    
    # ad transcript record
    new_records.append(transcript_record)
    
    # Add update exon records
    new_records.extend(exon_records_sorted)
    
    # add original CDS records (unchanged)
    new_records.extend(sorted(cds_records, key=lambda x: x.start))
    
    # add UTR records (only if they have positive length)
    if strand == '+':
        # 5' UTR
        if utr_5_start <= utr_5_end:
            utr_5_record = create_utr_record(first_exon, 'five_prime_utr', utr_5_start, utr_5_end)
            new_records.append(utr_5_record)
        
        # 3' UTR
        if utr_3_start <= utr_3_end:
            utr_3_record = create_utr_record(last_exon, 'three_prime_utr', utr_3_start, utr_3_end)
            new_records.append(utr_3_record)
    else:
        # 3' UTR (upstream for - strand)
        if utr_3_start <= utr_3_end:
            utr_3_record = create_utr_record(first_exon, 'three_prime_utr', utr_3_start, utr_3_end)
            new_records.append(utr_3_record)
        
        # 5' UTR (downstream for - strand)
        if utr_5_start <= utr_5_end:
            utr_5_record = create_utr_record(last_exon, 'five_prime_utr', utr_5_start, utr_5_end)
            new_records.append(utr_5_record)
    
    # Add any other records
    new_records.extend(other_records)
    
    return new_records

def main():
    parser = argparse.ArgumentParser(
        description="Add UTRs to GTF annotations where exon coordinates match CDS coordinates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gtf_add_utrs.py input.gtf output.gtf
  python gtf_add_utrs.py input.gtf output.gtf --utr5-length 200 --utr3-length 300
  python gtf_add_utrs.py input.gtf output.gtf -5 100 -3 150 --chrom-sizes chrom.sizes
        """
    )
    
    parser.add_argument('input_gtf', help='Input GTF file')
    parser.add_argument('output_gtf', help='Output GTF file')
    parser.add_argument('-5', '--utr5-length', type=int, default=150,
                        help='Length of 5\' UTR to add (default: 150)')
    parser.add_argument('-3', '--utr3-length', type=int, default=250,
                        help='Length of 3\' UTR to add (default: 250)')
    parser.add_argument('--chrom-sizes', type=str, default=None,
                        help='Chromosome sizes file (tab-separated: chromosome<tab>size)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.utr5_length < 0 or args.utr3_length < 0:
        print("Error: UTR lengths must be non-negative")
        sys.exit(1)
    
    input_file = args.input_gtf
    output_file = args.output_gtf
    utr_5_length = args.utr5_length
    utr_3_length = args.utr3_length
    chrom_sizes_file = args.chrom_sizes
    
    # Parse chromosome sizes if provided - will probably make mandatory in the future
    chrom_sizes = None
    if chrom_sizes_file:
        print(f"Reading chromosome sizes from: {chrom_sizes_file}")
        chrom_sizes = parse_chromosome_sizes(chrom_sizes_file)
        print(f"Loaded sizes for {len(chrom_sizes)} chromosomes")
    else:
        print("Warning: No chromosome sizes provided. UTRs may extend beyond chromosome boundaries.")
    
    print(f"Reading GTF file: {input_file}")
    print(f"UTR lengths: 5' = {utr_5_length}bp, 3' = {utr_3_length}bp")
    records = parse_gtf(input_file)
    
    print(f"Grouping {len(records)} records by transcript...")
    transcript_groups = group_by_transcript(records)
    
    print(f"Processing {len(transcript_groups)} transcripts...")
    
    processed_records = []
    modified_count = 0
    
    for transcript_id, transcript_records in transcript_groups.items():
        original_count = len(transcript_records)
        new_records = process_transcript(transcript_records, utr_5_length, utr_3_length, chrom_sizes)
        
        if len(new_records) > original_count:
            modified_count += 1
            print(f"  Added UTRs to transcript: {transcript_id}")
        
        processed_records.extend(new_records)
    
    # Add records that don't have transcript IDs (if any)
    for record in records:
        if not record.get_transcript_id():
            processed_records.append(record)
    
    print(f"Modified {modified_count} transcripts")
    print(f"Writing output to: {output_file}")
    
    # Sort records by chromosome and position
    processed_records.sort(key=lambda x: (x.seqname, x.start, x.end))
    
    with open(output_file, 'w') as f:
        for record in processed_records:
            f.write(record.to_line() + '\n')
    
    print("Done!")

if __name__ == "__main__":
    main()