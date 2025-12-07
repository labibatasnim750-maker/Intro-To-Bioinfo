
#!/bin/bash

# Change to the directory containing your FASTQ files
cd ~/Downloads/samples

# Create output directory for FASTA files
mkdir -p fasta_output

# Process each sample (A, B, C, D)
for sample in A B C D; do
    echo "Processing sample${sample}..."
    
    # Check if all three parts exist
    if [[ -f "sample${sample}_part1.FASTQ" && -f "sample${sample}_part2.FASTQ" && -f "sample${sample}_part3.FASTQ" ]]; then
        # Concatenate and convert FASTQ to FASTA
        # FASTQ format: 4 lines per record (@header, sequence, +, quality)
        # FASTA format: 2 lines per record (>header, sequence)
        
        # Method 1: Using awk (most efficient)
        cat "sample${sample}_part1.FASTQ" "sample${sample}_part2.FASTQ" "sample${sample}_part3.FASTQ" | \
        awk 'NR%4==1 {printf ">%s\n", substr($0,2)} NR%4==2 {print}' > "fasta_output/sample${sample}.fasta"
        
        # Alternative method using sed (also works)
        # cat "sample${sample}_part1.FASTQ" "sample${sample}_part2.FASTQ" "sample${sample}_part3.FASTQ" | \
        # sed -n '1~4s/^@/>/p;2~4p' > "fasta_output/sample${sample}.fasta"
        
        echo "Created: fasta_output/sample${sample}.fasta"
    else
        echo "Warning: Some files for sample${sample} are missing!"
    fi
done

echo "Done! All FASTA files are in the 'fasta_output' directory."


