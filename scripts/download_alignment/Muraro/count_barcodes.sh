#!/bin/bash
# Path to your FASTQ file and whitelist
fastq_file="Donor_D28_library_1_R1.fq.gz"
whitelist_file="/work/islet_cartography_scrna/whitelist/Muraro.whitelist"

# Count matches for the whitelist
# look at line 2 in fastq file
# extract character 1 to 8
revcount=$(zcat "$fastq_file" | awk 'NR % 4 == 2 {print substr($0, 1, 8)}' | grep -f "$whitelist_file" | wc -l)
echo "Matches for whitelist first 8 characters of fastq file 1 $whitelist_file: $revcount"

revcount=$(zcat "$fastq_file" | awk 'NR % 4 == 2 {print substr($0, 5, 8)}' | grep -f "$whitelist_file" | wc -l)
echo "Matches for whitelist for characters on position 5 to 18 of fastq file 1 $whitelist_file: $revcount"
