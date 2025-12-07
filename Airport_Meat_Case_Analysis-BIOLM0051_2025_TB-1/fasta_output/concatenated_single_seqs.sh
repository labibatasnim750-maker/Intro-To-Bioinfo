#!/bin/bash

# Process each sample
for sample in sampleA sampleB sampleC sampleD; do
    input_file="${sample}_cleaned.fasta"
    output_file="${sample}.fasta"
    
    if [ -f "$input_file" ]; then
        echo "Processing $input_file..."
        
        # Count parts
        parts=$(grep -c ">" "$input_file")
        echo "  Found $parts parts in $input_file"
        
        # Combine into one sequence
        grep -v ">" "$input_file" | tr -d '\n' > temp_seq.txt
        
        # Create new file with single header
        echo ">$sample" > "$output_file"
        cat temp_seq.txt >> "$output_file"
        echo "" >> "$output_file"
        
        # Check length
        length=$(grep -v ">" "$output_file" | wc -c)
        echo "  Combined length: $length bases"
        
        rm temp_seq.txt
    else
        echo "Warning: $input_file not found"
    fi
done

echo "Done! Created single-sequence files:"
ls -la sample?.fasta 2>/dev/null || ls sampleA.fasta sampleB.fasta sampleC.fasta sampleD.fasta 2>/dev/null
