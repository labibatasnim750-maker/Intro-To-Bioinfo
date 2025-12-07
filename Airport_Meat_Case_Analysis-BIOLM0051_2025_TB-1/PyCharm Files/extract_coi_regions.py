#!/usr/bin/env python3
"""
EXTRACTING COI GENE REGION FROM REFERENCE GENOMES
Matching the region corresponding to the sample sequences
"""
# Importing...
from pathlib import Path
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

print(" EXTRACTING COI GENE REGIONS FROM REFERENCE GENOMES")
print("=" * 60)

# Paths
home = Path.home()
downloads = home / "Downloads"
samples_dir = downloads / "samples" / "fasta_output"
refs_dir = downloads

# Creating output directory
output_dir = Path.cwd() / "extracted_coi_regions"
output_dir.mkdir(exist_ok=True)
print(f" Output directory: {output_dir}")


def find_coi_region(sample_seq, ref_seq, sample_id=""):
    """
    Finding the COI gene region in reference sequence that matches the sample.
    """
    sample_len = len(sample_seq)

    # COI gene specific patterns (conserved in vertebrates)
    coi_patterns = [
        "ATGTTC",  # Common COI start in vertebrates
        "CATAA",  # Another common pattern
        "GGGGCA",  # From your sampleC
        "CTCGCT",  # From your sampleD
        "TTGGA",  # Common COI pattern
        "CCTGG",  # Common COI pattern
        "GGAGG",  # Common COI pattern
        "CAAAC",  # Common COI pattern
    ]

    best_match = None
    best_score = 0
    best_pattern = ""

    # Searching for patterns in the reference
    for pattern in coi_patterns:
        pos = ref_seq.find(pattern)
        if pos != -1:
            # Extracting region similar to sample length
            start = max(0, pos - 50)  # Some buffer before pattern
            end = min(len(ref_seq), start + sample_len + 100)  # Sample length + buffer

            extracted = ref_seq[start:end]

            # Simple scoring: count matching bases in first 100 positions
            compare_len = min(100, len(sample_seq), len(extracted))
            score = 0
            for i in range(compare_len):
                if sample_seq[i] == extracted[i]:
                    score += 1

            # Bonus for longer matches
            if score > best_score:
                best_score = score
                best_match = (start, end, extracted)
                best_pattern = pattern

    return best_match, best_pattern, best_score


def main():
    all_extracted_records = []
    all_combined_records = []

    for sample_letter in ["A", "B", "C", "D"]:
        print(f"\n PROCESSING SAMPLE {sample_letter}")
        print("-" * 40)

        # Reading sample sequence
        sample_file = samples_dir / f"sample{sample_letter}.fasta"
        if not sample_file.exists():
            print(f" Sample file not found: {sample_file}")
            continue

        sample_records = list(SeqIO.parse(sample_file, "fasta"))
        if not sample_records:
            print(f" No sequences in {sample_file}")
            continue

        sample_record = sample_records[0]
        sample_seq = str(sample_record.seq).upper()
        print(f"Sample: {sample_record.id}")
        print(f"Length: {len(sample_seq)} bp")
        print(f"First 30 bases: {sample_seq[:30]}")

        # Reading reference sequences
        ref_file = refs_dir / f"refseq_for_sample{sample_letter}.txt"
        if not ref_file.exists():
            print(f" Reference file not found: {ref_file}")
            continue

        ref_records = list(SeqIO.parse(ref_file, "fasta"))
        print(f"Found {len(ref_records)} reference sequences")

        extracted_records = []

        for ref_record in ref_records:
            ref_seq = str(ref_record.seq).upper()
            ref_len = len(ref_seq)

            print(f"\n  Reference: {ref_record.id}")
            print(f"    Original length: {ref_len:,} bp")

            # Trying to extract COI region
            result, pattern, score = find_coi_region(sample_seq, ref_seq, ref_record.id)

            if result and score > 20:  # Reasonable match threshold
                start, end, extracted_seq = result
                print(f"    ✓ Found pattern '{pattern}' at position {start}")
                print(f"    Extracted region: {start:,}-{end:,} ({len(extracted_seq)} bp)")
                print(f"    Match score: {score}/100")

                # Creating new record
                new_id = f"{ref_record.id}_COI_extract"
                new_desc = f"COI region from {ref_record.description} | pos:{start}-{end} | score:{score}"

                new_record = SeqRecord(
                    Seq.Seq(extracted_seq),
                    id=new_id,
                    description=new_desc
                )

                extracted_records.append(new_record)
                all_extracted_records.append(new_record)

            else:
                print(f"    Could not find good COI match (score: {score})")
                print(f"    Using first {len(sample_seq)} bp as fallback")

                # Using beginning of sequence
                extracted_seq = ref_seq[:len(sample_seq) + 100]
                new_id = f"{ref_record.id}_first{len(extracted_seq)}bp"
                new_desc = f"First {len(extracted_seq)} bp of {ref_record.id} (fallback)"

                new_record = SeqRecord(
                    Seq.Seq(extracted_seq),
                    id=new_id,
                    description=new_desc
                )

                extracted_records.append(new_record)
                all_extracted_records.append(new_record)

        # Saving extracted references for this sample
        if extracted_records:
            output_file = output_dir / f"sample{sample_letter}_extracted_refs.fasta"
            SeqIO.write(extracted_records, output_file, "fasta")
            print(f"\n  Saved {len(extracted_records)} extracted references to:")
            print(f"     {output_file.name}")

            # Saving with samples for combined file
            combined_file = output_dir / f"sample{sample_letter}_with_refs.fasta"
            combined_records = [sample_record] + extracted_records
            SeqIO.write(combined_records, combined_file, "fasta")
            print(f"   Saved combined sample + references to:")
            print(f"     {combined_file.name}")

            # Adding to overall combined records
            all_combined_records.extend(combined_records)

    # Saving all extracted sequences together
    if all_extracted_records:
        all_extracted_file = output_dir / "all_extracted_references.fasta"
        SeqIO.write(all_extracted_records, all_extracted_file, "fasta")
        print(f"\n All extracted references saved to: {all_extracted_file.name}")

        # Showing statistics
        lengths = [len(rec.seq) for rec in all_extracted_records]
        print(f" Statistics for extracted references:")
        print(f"  Total sequences: {len(all_extracted_records)}")
        print(f"  Average length: {sum(lengths) / len(lengths):.1f} bp")
        print(f"  Min length: {min(lengths)} bp")
        print(f"  Max length: {max(lengths)} bp")

    # Saving complete combined file
    if all_combined_records:
        complete_file = output_dir / "complete_dataset.fasta"
        SeqIO.write(all_combined_records, complete_file, "fasta")
        print(f"\n Complete dataset (samples + references) saved to: {complete_file.name}")

        # Counting samples vs references
        sample_count = sum(1 for rec in all_combined_records if 'sample' in rec.id.lower())
        ref_count = len(all_combined_records) - sample_count
        print(f"  • Sample sequences: {sample_count}")
        print(f"  • Reference sequences: {ref_count}")

    print("\n" + "=" * 60)
    print(" NEXT STEPS:")
    print("=" * 60)
    print("1. Check the extracted sequences in 'extracted_coi_regions' folder")
    print("2. Verify sequences are ~1,400-1,600 bp (similar to your samples)")
    print("3. Translate these extracted sequences")
    print("4. Align with MUSCLE")
    print(f"\n Your main file for translation: {output_dir}/complete_dataset.fasta")


if __name__ == "__main__":
    main()