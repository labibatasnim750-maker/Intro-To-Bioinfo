#!/usr/bin/env python3
"""
FINAL COI TRANSLATION SCRIPT
"""
# Importing...
from pathlib import Path
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

print(" FINAL COI TRANSLATION")
print("=" * 60)

# Working Directory: /Users/labibatasnim/PycharmProjects/HelloWorld/extracted_coi_regions/
base_dir = Path("/Users/labibatasnim/PycharmProjects/HelloWorld")
extracted_dir = base_dir / "extracted_coi_regions"
input_file = extracted_dir / "complete_dataset.fasta"
output_dir = base_dir / "final_coi_proteins"
output_dir.mkdir(exist_ok=True)

print(f" Looking for: {input_file}")

if not input_file.exists():
    print(f" ERROR: Input file not found!")
    print(f"\nAvailable files in {extracted_dir}:")
    for file in extracted_dir.glob("*"):
        if file.is_file():
            print(f"  • {file.name}")
    exit()

print(f" Found input file!")
print(f" Output directory: {output_dir}")

# Reading sequences
records = list(SeqIO.parse(input_file, "fasta"))
print(f" Total sequences: {len(records)}")

# Showing breakdowns
sample_seqs = [rec for rec in records if 'sample' in rec.id.lower()]
ref_seqs = [rec for rec in records if 'COI' in rec.id or 'extract' in rec.id.lower()]
other_seqs = len(records) - len(sample_seqs) - len(ref_seqs)

print(f" Breakdown:")
print(f"  • Sample sequences: {len(sample_seqs)}")
print(f"  • Reference COI sequences: {len(ref_seqs)}")
if other_seqs > 0:
    print(f"  • Other sequences: {other_seqs}")

print(f"\n Using Genetic Code 2 (Vertebrate Mitochondrial)")
genetic_code = 2


def find_best_protein(dna_seq):
    """
    Finding the best protein translation by trying all 6 reading frames.
    """
    dna_seq_obj = Seq.Seq(dna_seq)
    best_protein = None
    best_length = 0
    best_frame = 0
    best_strand = "forward"

    # Trying forward strand (3 frames)
    for frame in [0, 1, 2]:
        try:
            # Extracting frame
            frame_seq = dna_seq_obj[frame:]
            # Making multiple of 3
            if len(frame_seq) % 3 != 0:
                frame_seq = frame_seq[:-(len(frame_seq) % 3)]

            if len(frame_seq) < 100:
                continue

            # Translating...
            protein_seq = frame_seq.translate(table=genetic_code, to_stop=True)

            if len(protein_seq) > best_length:
                best_protein = protein_seq
                best_length = len(protein_seq)
                best_frame = frame
                best_strand = "forward"
        except:
            continue

    # Trying reverse complement (3 frames)
    rev_seq = dna_seq_obj.reverse_complement()
    for frame in [0, 1, 2]:
        try:
            frame_seq = rev_seq[frame:]
            if len(frame_seq) % 3 != 0:
                frame_seq = frame_seq[:-(len(frame_seq) % 3)]

            if len(frame_seq) < 100:
                continue

            protein_seq = frame_seq.translate(table=genetic_code, to_stop=True)

            if len(protein_seq) > best_length:
                best_protein = protein_seq
                best_length = len(protein_seq)
                best_frame = frame
                best_strand = "reverse"
        except:
            continue

    return best_protein, best_length, best_frame, best_strand


print(f"\n TRANSLATING SEQUENCES...")
print("-" * 60)

translated_records = []
stats = []

for i, record in enumerate(records, 1):
    dna_seq = str(record.seq)
    dna_len = len(dna_seq)

    print(f"\n{i:2d}. {record.id}")
    print(f"    DNA: {dna_len} bp")

    # Finding best protein
    protein, protein_len, frame, strand = find_best_protein(dna_seq)

    if protein and protein_len >= 100:  # Minimum 100 aa
        print(f"     Protein: {protein_len} aa")
        print(f"       Frame: {frame}, Strand: {strand}")

        # Checking if reasonable COI length
        expected_aa = dna_len // 3
        if abs(protein_len - expected_aa) > 100:
            print(f"         Length discrepancy: expected ~{expected_aa} aa")

        # Creating new record
        new_id = f"{record.id}_COI_protein"
        new_desc = f"COI protein | DNA:{dna_len}bp | AA:{protein_len} | frame:{frame} | strand:{strand}"

        new_record = SeqRecord(
            protein,
            id=new_id,
            description=new_desc
        )

        translated_records.append(new_record)
        stats.append((record.id, dna_len, protein_len, strand))

    else:
        print(f"     No good protein found (best: {protein_len if protein else 0} aa)")

# Saving results
if translated_records:
    output_file = output_dir / "coi_proteins_aligned.fasta"
    SeqIO.write(translated_records, output_file, "fasta")

    print(f"\n" + "=" * 60)
    print(f" TRANSLATION COMPLETE!")
    print("=" * 60)

    print(f"\n Saved to: {output_file}")
    print(f" Successfully translated: {len(translated_records)}/{len(records)} sequences")

    # Statistics
    protein_lengths = [s[2] for s in stats]
    print(f"\n PROTEIN LENGTH STATISTICS:")
    print(f"  Minimum: {min(protein_lengths)} aa")
    print(f"  Maximum: {max(protein_lengths)} aa")
    print(f"  Average: {sum(protein_lengths) / len(protein_lengths):.1f} aa")

    # Counting by strand
    forward_count = sum(1 for s in stats if s[3] == "forward")
    reverse_count = len(stats) - forward_count
    print(f"\n STRAND DISTRIBUTION:")
    print(f"  Forward strand: {forward_count}")
    print(f"  Reverse strand: {reverse_count}")

    # Showing results table
    print(f"\n RESULTS:")
    print("-" * 70)
    print(f"{'Sequence':<30} {'DNA':>6} {'Protein':>8} {'Strand':>8}")
    print("-" * 70)

    for seq_id, dna_len, protein_len, strand in stats:
        short_id = seq_id[:27] + "..." if len(seq_id) > 30 else seq_id
        print(f"{short_id:<30} {dna_len:>6} {protein_len:>8} {strand:>8}")

    print(f"\n NEXT STEP: ALIGN WITH MUSCLE")
    print("=" * 60)
    print(f"\nFile to align: {output_file}")

    # Checking for stop codons
    print(f"\n  QUALITY CHECK:")
    internal_stops = 0
    for rec in translated_records:
        if '*' in str(rec.seq):
            internal_stops += 1

    if internal_stops == 0:
        print(f"   No sequences have stop codons")
    else:
        print(f"    {internal_stops} sequences have stop codons")

else:
    print(f"\n NO PROTEINS TRANSLATED!")
    print(f"\nTrying alternative approach...")

    # Trying a different method
    print(f"\n ALTERNATIVE: Try extracting proteins from beginning")

    simple_translated = []
    for record in records:
        dna_seq = Seq.Seq(str(record.seq))

        # Trying to translate from beginning
        try:
            # Making multiple of 3
            remainder = len(dna_seq) % 3
            if remainder != 0:
                dna_seq = dna_seq[:-remainder]

            protein = dna_seq.translate(table=genetic_code, to_stop=True)

            if len(protein) > 100:
                print(f"  ✓ {record.id}: {len(protein)} aa")
                new_record = SeqRecord(
                    protein,
                    id=f"{record.id}_simple",
                    description=f"Simple translation | {len(protein)} aa"
                )
                simple_translated.append(new_record)
        except:
            pass

    if simple_translated:
        output_file = output_dir / "coi_proteins_simple.fasta"
        SeqIO.write(simple_translated, output_file, "fasta")
        print(f"\n Saved {len(simple_translated)} proteins to: {output_file}")
    else:
        print(f"\n Even simple translation failed!")