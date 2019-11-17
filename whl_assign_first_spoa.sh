# make bedfile from cmh file: ~/RESEARCH/SpOA/combined_data/analysis/cmh.out

cat ~/RESEARCH/SpOA/combined_data/analysis/cmh.out.txt | tail -n +2 | awk '{OFS="\t"} { print $1, $2-1, $2 }' | sort -k1,1 -k2,2n > ~/RESEARCH/SpOA/combined_data/analysis/cmh.all.sorted.bed

cat ~/RESEARCH/SpOA/combined_data/analysis/cmh.out.txt  |sort -k1,1 -k2,2n > ~/RESEARCH/SpOA/combined_data/analysis/cmh.master.sort.out


### assign whl names to genes

# all cmh
/users/r/b/rbrennan/bin/bedtools2/bin/bedtools closest -a ~/RESEARCH/SpOA/combined_data/analysis/cmh.all.sorted.bed \
-b /users/r/b/rbrennan/reference/urchin_probes.sorted.bed \
-d -t first  > ~/RESEARCH/SpOA/combined_data/analysis/cmh.all.genes.txt

# convert cmh.genes.txt gene names to those compatible with echinobase
#get id from 5th column

# remove transcript number from whl names
cat ~/RESEARCH/SpOA/combined_data/analysis/cmh.all.genes.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
    cut -f 1-2 -d '.' > ~/RESEARCH/SpOA/combined_data/analysis/cmh.all.id
