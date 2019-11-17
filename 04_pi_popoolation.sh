#!/bin/bash

#sort the bam file
# D1_control.bam is D1 pH 8.1S
# Other samples: D7_control.bam, D7_7_5S.bam, D7_7_0S.bam, D7_7_5V.bam, D7_7_0V.bam, D7_8_1V.bam

# Need to merge bam files accordingly first & sort (done)

#/users/r/b/rbrennan/bin/samtools-1.6/samtools view -q 20 -b ~/RESEARCH/SpOA/combined_data/mapped/D1.bam | /users/r/b/rbrennan/bin/samtools-1.6/samtools sort - -o ~/RESEARCH/SpOA/combined_data/mapped/D1.sort.bam
#/users/r/b/rbrennan/bin/samtools-1.6/samtools view -q 20 -b ~/RESEARCH/SpOA/combined_data/mapped/D7_8.bam | /users/r/b/rbrennan/bin/samtools-1.6/samtools sort -  -o ~/RESEARCH/SpOA/combined_data/mapped/D7_8.sort.bam
#/users/r/b/rbrennan/bin/samtools-1.6/samtools view -q 20 -b ~/RESEARCH/SpOA/combined_data/mapped/D7_7.bam | /users/r/b/rbrennan/bin/samtools-1.6/samtools sort - -o ~/RESEARCH/SpOA/combined_data/mapped/D7_7.sort.bam

#convert to pileup:
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D1_control.bam > ~/RESEARCH/SpOA/combined_data/mapped/D1_control.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_control.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_control.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5V.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5V.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0V.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0V.pileup
/users/r/b/rbrennan/bin/samtools-1.6/samtools mpileup ~/RESEARCH/SpOA/combined_data/mapped/D7_8_1V.bam > ~/RESEARCH/SpOA/combined_data/mapped/D7_8_1V.pileup

# calc pi
perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D1_control.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D1_control.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20  --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_control.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_control.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20  --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5V.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5V.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0V.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0V.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_8_1V.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_8_1V.pi --measure pi --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000


# tajima's D

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D1_control.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D1_control.TD --measure D --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_5S.TD --measure D --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000

perl /users/r/b/rbrennan/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.pileup --output ~/RESEARCH/SpOA/combined_data/mapped/D7_7_0S.TD --measure D --window-size 400 --step-size 200 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 20 --min-covered-fraction 0.5 --pool-size 1000