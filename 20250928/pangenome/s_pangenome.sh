# Rebecca Satterwhite
# 4/11/25 anvio on 'big' pangenome for 5 pathovar manuscript

############################################################
# # # # # # # # # install
conda config --env --set subdir osx-64
conda activate /Users/.../opt/anaconda3/envs/anvio-8

conda install -y -c conda-forge -c bioconda python=3.10 \
        sqlite=3.46 prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript \
        nodejs=20.12.2

conda install -y -c bioconda fastani

curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz

export CC=/usr/bin/clang
export CXX=/usr/bin/clang++

pip install anvio-8.tar.gz # all good
anvi-self-test --suite mini # all good

anvi-setup-scg-taxonomy
anvi-setup-ncbi-cogs
anvi-setup-kegg-data
anvi-self-test --suite pangenomics

# got a PyANI error, so ran the following 
pip install matplotlib==3.5.1
anvi-self-test --suite pangenomics # all god

# skipped 6.2 installation bc not using metagenomes

# needed this later (4/14)
conda install llvmlite
############################################################

# # # clean & QC genomes, build pangenome
conda activate /.../opt/anaconda3/envs/anvio-8

for sample in `cat list`; do anvi-script-reformat-fasta $sample -o $sample.fix  -l 0 --simplify-names; done
# manually rename to get rid of the original extension

for sample in *fna; do anvi-script-reformat-fasta $sample -o $sample.fix -l 0 --simplify-names; done
anvi-script-reformat-fasta CFBP3846.fna -o CFBP3846.fix -l 0 --simplify-names --seq-type NT
anvi-script-reformat-fasta ICMP14479.fna -o ICMP14479.fix -l 0 --simplify-names


# Heads up that this takes forever
for sample in *.fix; do anvi-gen-contigs-database -f $sample -o $sample.db; done
anvi-gen-contigs-database -f CFBP3846.fix -o CFBP3846.fix.db


# call single copy core genes - must do this to estimate completeness
for sample in *db; do anvi-run-hmms -c $sample; done


# summary stats
anvi-estimate-genome-completeness -e external-genomes.txt -o completition.txt


# functional annotation 
anvi-setup-ncbi-cogs --num-threads 10 --cog-data-dir /Users/.../big_pg/C20 --cog-version COG20
#anvi-setup-pfams --pfam-data-dir /home/rsatterwhite/ 

# Heads up that this takes forever
for sample in *db; do anvi-run-ncbi-cogs -c $sample --cog-data-dir /Users/.../big_pg/C20 --cog-version COG20 --num-threads 10; done
#for sample in *db; do anvi-run-pfams -c $sample --pfam-data-dir /home/rsatterwhite/; done


# build a genomes storage
# names for external-genomes dont matter/can't contain . or start with a # (put As in front of #s, removed _)
anvi-gen-genomes-storage -e external-genomes.txt -o big37-GENOMES.db


# mcl inflation parameter because we are looking at strains of the same species/close relatives
anvi-pan-genome -g big37-GENOMES.db --project-name "big37_pan" --num-threads 10 --minbit 0.5 --mcl-inflation 10

anvi-display-pan -p big37_pan/big37_pan-PAN.db -g big37-GENOMES.db

cd /Users/.../2025/pg_small
anvi-display-pan -p small5_pan/small5_pan-PAN.db -g small5-GENOMES.db

############################################################
# operations/additional analyses
conda activate /Users/.../opt/anaconda3/envs/anvio-8
cd /Users/.../2025/big_pg


#ANI
anvi-compute-genome-similarity --external-genomes external-genomes.txt \
--program fastANI \
--output-dir fastANI \
--num-threads 20 \
--pan-db big37_pan/big37_pan-PAN.db


# extract singleton genes (DNA seqs (flag) or default is protein)
anvi-get-sequences-for-gene-clusters -g big37-GENOMES.db \
                                     -p big37_pan/big37_pan-PAN.db \
                                     -o singleton-gene-prot-fasta \
                                     --max-num-genes-from-each-genome 1 \
                                     --max-num-genomes-gene-cluster-occurs 1 
                                     
                                     --report-DNA-sequences --force-overwrite
# # # # deviating from anvio here
# pull only the singeltons for the focal 5
awk '/^>/ {f=($0 ~ /A1448A/)} f' singleton-gene-nt-2.fasta > singletons_1448A.fasta
awk '/^>/ {f=($0 ~ /A9/)} f' singleton-gene-nt-2.fasta > singletons_A9.fasta
awk '/^>/ {f=($0 ~ /DC3000/)} f' singleton-gene-nt-2.fasta > singletons_DC.fasta
awk '/^>/ {f=($0 ~ /ES4326/)} f' singleton-gene-nt-2.fasta > singletons_ES.fasta
awk '/^>/ {f=($0 ~ /NP29/)} f' singleton-gene-nt-2.fasta > singletons_NP.fasta

# expected singleton genes (not exact but close approximation)
178 a14 / 71 a9 / 125 dc / 671 es / 66 np

# seqs in a fasta
grep -c "^>" singletons_1448A.fasta # 175
grep -c "^>" singletons_A9.fasta # 71
grep -c "^>" singletons_DC.fasta # 121
grep -c "^>" singletons_ES.fasta # 658
grep -c "^>" singletons_NP.fasta # 66


# pull only the singeltons for the focal 5
awk '/^>/ {f=($0 ~ /A1448A/)} f' singleton-gene-prot.fasta > singletons_prot_1448A.fasta
awk '/^>/ {f=($0 ~ /A9/)} f' singleton-gene-prot.fasta > singletons_prot_A9.fasta
awk '/^>/ {f=($0 ~ /DC3000/)} f' singleton-gene-prot.fasta > singletons_prot_DC.fasta
awk '/^>/ {f=($0 ~ /ES4326/)} f' singleton-gene-prot.fasta > singletons_prot_ES.fasta
awk '/^>/ {f=($0 ~ /NP29/)} f' singleton-gene-prot.fasta > singletons_prot_NP.fasta


# run bakta on singleton fastas on galaxy, download annotations + tabular summary file
cd /Users/.../2025/singletons

# run custom py script to add the COG annotation for each seq
python script.py path/to/input.tabular path/to/cog_database.tsv path/to/output.tabular

# sadly there aren't many COG categories called
python s_clean-anns.py sings-A14.tabular clean-sings-A14.tabular
python s_clean-anns.py sings-ES.tabular cog_database.tsv clean-sings-ES.tabular
python s_clean-anns.py sings-DC.tabular cog_database.tsv clean-sings-DC.tabular
python s_clean-anns.py sings-A9.tabular cog_database.tsv clean-sings-A9.tabular
python s_clean-anns.py sings-NP.tabular cog_database.tsv sings-clean-NP.tabular
