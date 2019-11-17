# Script for combining the merged.bam files
```
/data/programs/samtools merge complete.merged.bam *.bam
```
Screen ID: 7159.pts-0.pespenilab

# Making a script for Sorting & Indexing the complete merged bam file

# SortIndexSamtools.sh in /users/a/m/amakukho/RESEARCH/SpOA/scripts
```
#!/bin/bash

cd /users/a/m/amakukho/RESEARCH/SpOA/combined_data/mapped

## Sorting
# NOTE: Need to specify a different tmp directory or it will fill the rhel-root...
/data/programs/sambamba_v0.6.0 sort -m 8G -t 8 -p --tmpdir=/users/a/m/amakukho/RESEARCH/SpOA/combined_data/tmp complete.merged.bam complete.merged.sorted.bam

## Indexing
/data/programs/samtools index complete.merged.sorted.bam
```

##??? At the end of the screen: [W::bam_hdr_read] EOF marker is absent. The input is probably truncated.

--------------------------------------------------------------------------------------------------------
• When I was using the script to count reads in the bamfile, I thought it was counting only the mapped reads, but it turns out bam files have both the mapped and unmapped reads in them.

Reid had suggested using samtools flags to determine how many reads were mapped, can do this with:

```
samtools view -c -F 4 BAMFILE

```

# This code below will count the high quality mapped reads
```
## Example
samtools view -c -F 4 -q 20 bam

## Code used
/data/programs/samtools view -c -F 4 -q 20 complete.merged.sorted.bam
```


# When running the above code, got the following message, going to look into what this means...
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated.

[main_samview] truncated file.


# Tried running the same code but on the only merged bam file (before being sorted with sambamba), and I did not get the same error as above.

```
/data/programs/samtools view -c -F 4 -q 20 complete.merged.bam
```




# Some insight from Seq Answers:
The cryptic error from samtools "EOF marker is absent" is referring to the absence of a special empty BGZF block of 28 bytes, which samtools looks for at the end of the data to indicate the BAM file is complete.

If you see that error, either:

(a) Your file is somehow truncated or incomplete (a real error)
(b) Your file is from a tool not writing this EOF marker (perhaps a very old samtools?)

# To address the above issue, we can look for the EOF marker at the bottom that we'd expect for the final 28 bytes of a completed bam file using the following code:

```
tail BAMFILE | hexdump -C
```
Starts with 1f, End in 9 sets of 00

# The bam files before sorting don't have this issue, so it appears sorting is what caused this issue. Now going to check the total number of mapped reads in the complete.merged.bam file (pre-sorting)
```
/data/programs/samtools view -c -F 4 -q 20 complete.merged.bam
```


--------------------------------------------------------------------------------------------------------

                            ###  CHECKING ON MAPPING STATUS OF FILES ###

# Ran the following command to look at percent mapping for the complete.merged.bam file
```
/data/programs/samtools flagstat complete.merged.bam
```
Output:
958045418 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
127217774 + 0 supplementary
0 + 0 duplicates
956006630 + 0 mapped (99.79%:-nan%)
830827644 + 0 paired in sequencing
415413822 + 0 read1
415413822 + 0 read2
588332692 + 0 properly paired (70.81%:-nan%)
826750068 + 0 with itself and mate mapped
2038788 + 0 singletons (0.25%:-nan%)
230124044 + 0 with mate mapped to a different chr
122933809 + 0 with mate mapped to a different chr (mapQ>=5)

# Interpretation

956,006,630 reads mapped total (which was 99.79% of the reads). Of those, 614,297,893 reads had a Q20 or higher score (good quality). This means ~64.26% of the mapped reads were high quality


--------------------------------------------------------------------------------------------------------

## Indexing
```
/data/programs/samtools index complete.merged.sorted.bam
```
