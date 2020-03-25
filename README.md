# kAAmer analyses with demonstration

This project includes python scripts for different bacterial genomics / metagenomics analyses.

You can follow this guide as a demo for typical bacterial genomics analyses which uses kAAmer has the database engine.

You first need to install kAAmer: see https://zorino.github.io/kaamer.

You should also use conda / miniconda for the python dependencies : see https://docs.conda.io/en/latest/miniconda.html.

Another project, ([microbe-dbs](https://github.com/zorino/microbe-dbs.git)), is used to download relevant databases.

``` shell
# Clone this repo
git clone https://github.com/zorino/kaamer_analyses.git

# Create conda environment
conda env create -f  kaamerpy-env.yml
conda activate kaamerpy
```


## Antibiotic resistance gene identification

For a demonstration purpose we uploaded in this repo a pan-resistant strain of Pseudomonas
aeruginosa [E6130952](https://www.ncbi.nlm.nih.gov/biosample/SAMN06349407) ([chromosome](https://www.ncbi.nlm.nih.gov/nuccore/CP020603.1) + [plasmid](https://www.ncbi.nlm.nih.gov/nuccore/CP020602.1)).

This example download the latest ncbi ARG database, build it using kAAmer and identify resistance
genes in the E6130952 strain.


``` shell
# Go to this repo directory
cd kaamer_analyses

# Download ncbi_arg database..
git clone https://github.com/zorino/microbe-dbs.git

## note that you will need python and the bio lib to run this:
microbe-dbs/microbe-dbs ncbi_arg demo-data/ncbi_arg

# Create the kAAmer database
kaamer-db -make -f tsv -i demo-data/ncbi_arg/ncbi-arg_$(date +%Y-%m-%d)/ReferenceGeneCatalog.tsv -d demo-data/ncbi_arg.kaamer

# Search for ARG in genome with kaamer
kaamer-db -server -d demo-data/ncbi_arg.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t nt -i demo-data/Pae_E6130952.fa -o demo-data/Pae_E6130952.arg.tsv -aln -ann

# Analyse ARG Results
python arg_identifier.py demo-data/Pae_E6130952.arg.tsv > demo-data/Pae_E6130952.arg.sum.tsv

# Check the ARG results: demo-data/Pae_E6130952.arg.sum.tsv
less demo-data/Pae_E6130952.arg.sum.tsv

```

## Bacterial genome annotation

We are going to use the chromosome of the Pseudomonas aeruginosa E6130952 strain (CP020603.1) for automatic
annotation using kAAmer.

``` shell
# Go to this repo directory
cd kaamer_analyses

# Download prebuilt kAAmer database for the Pseudomonadaceae family.
wget https://kaamer.genome.ulaval.ca/dwl/Pseudomonadaceae.kaamer.tgz -O demo-data/Pseudomonadaceae.kaamer.tgz
tar xvf demo-data/Pseudomonadaceae.kaamer.tgz -C demo-data/

# Search for ARG in genome with kaamer
kaamer-db -server -d demo-data/Pseudomonadaceae.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t nt -i demo-data/CP020603.1.fa -o demo-data/CP020603.1.fa.ann.tsv -aln -ann

# Create GFF annotation file
python genome_annotation.py --seq demo-data/CP020603.1.fa --kaamer_res demo-data/CP020603.1.fa.ann.tsv > demo-data/CP020603.1.gff

```


## Metagenomics profiling (gut metagenome)

For the metagenomic demo, we are going to download one read from a gut metagenome publish by Nielsen
et al.[1] and make the profiling based on the Mgnify database of the human gut [2].

Downloading the uhgp-90 database can take a while so expect this demo to be longer than usual.

Also running this kAAmer database require a fair amount of RAM, plan for at least 16 GB.

``` shell
# Go to this repo directory
cd kaamer_analyses

# Download prebuilt kAAmer database of Mgnify for gut metagenome profiling
wget https://kaamer.genome.ulaval.ca/dwl/uhgp-90_2019-09.kaamer.tgz -O demo-data/uhgp-90_2019-09.kaamer.tgz
tar xvf demo-data/uhgp-90-2019_09.kaamer.tgz -C demo-data/

# Download demo gut metagenome fastq 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR209/ERR209293/ERR209293_1.fastq.gz -O demo-data/mg-reads.fastq.gz

# Start the Mgnify Human Gut kaamer database...
kaamer-db -server -d demo-data/uhgp-90_2019-09.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t fastq -m 3 -i demo-data/mg-reads.fastq.gz -o demo-data/mg-reads.ann.tsv -ann

# Analyse 
python metagenome_profiling.py demo-data/mg-reads.fastq.gz demo-data/mg-reads.ann.tsv demo-data/mg-reads-results

```

[1] Nielsen, H. B., Almeida, M., Juncker, A. S., Rasmussen, S., Li, J., Sunagawa, S., … MetaHIT Consortium. (2014). Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes. Nature Biotechnology, 32(8), 822–828. https://doi.org/10.1038/nbt.2939

[2] Mitchell, A. L., Almeida, A., Beracochea, M., Boland, M., Burgin, J., Cochrane, G., … Finn, R. D. (2019). MGnify: the microbiome analysis resource in 2020. Nucleic Acids Research. https://doi.org/10.1093/nar/gkz1035
