# kAAmer bacterial analyses

Bacterial genomics analyses demonstration which uses kAAmer has the database engine.

You first need to install kAAmer: see https://zorino.github.io/kaamer.

Another project, ([microbe-dbs](https://github.com/zorino/microbe-dbs.git)), is used to download relevant databases.

``` shell
# Clone this repo
git clone https://github.com/zorino/kaamer-demo.git
cd kaamer-bacterial-analyses
```

## Antibiotic Resistance Gene identification

For a demonstration purpose we uploaded in this repo a pan-resistant strain of Pseudomonas
aeruginosa [E6130952](https://www.ncbi.nlm.nih.gov/biosample/SAMN06349407) ([chromosome](https://www.ncbi.nlm.nih.gov/nuccore/CP020603.1) + [plasmid](https://www.ncbi.nlm.nih.gov/nuccore/CP020602.1)).

This example download the latest ncbi ARG database, build it using kAAmer and identify resistance
genes in the E6130952 strain.


``` shell
# Go to this repo directory
cd kaamer-demo

# Download ncbi_arg database..
git clone https://github.com/zorino/microbe-dbs.git

## note that you will need python and the bio lib to run this:
microbe-dbs/microbe-dbs ncbi_arg data/ncbi_arg

# Create the kAAmer database
kaamer-db -make -f tsv -i data/ncbi_arg/ncbi-arg_$(date +%Y-%m-%d)/ReferenceGeneCatalog.tsv -d data/ncbi_arg.kaamer

# Search for ARG in genome with kaamer
kaamer-db -server -d data/ncbi_arg.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t nt -i data/Pae_E6130952.fsa -o data/Pae_E6130952.arg.tsv -aln -ann

# Analyse ARG Results
python scripts/arg-identifier.py data/Pae_E6130952.arg.tsv > data/Pae_E6130952.arg.sum.tsv

```

## Bacterial Genome Annotation

We are going to use the chromosome of the Pseudomonas aeruginosa E6130952 strain (CP020603.1) for automatic
annotation using kAAmer.

``` shell
# Go to this repo directory
cd kaamer-demo

# Download prebuilt kAAmer database for the Pseudomonadaceae family.
wget xxx

# Search for ARG in genome with kaamer
kaamer-db -server -d data/Pseudomonadaceae.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t nt -i data/CP020603.1.fa -o data/CP020603.1.fa.ann.tsv -aln -ann

# Create GFF annotation file
python scripts/genome-annotation.py --seq data/CP020603.1.fa --kaamer_res data/CP020603.1.fa.ann.tsv > data/CP020603.1.gff

```


## Metagenomics profiling (gut metagenome)

For the metagenomic demo, we are going to download one read from a gut metagenome publish by Nielsen
et al.[^1] and make the profiling based on the Mgnify database of the human gut [^2].

[^1]: Nielsen, H. B., Almeida, M., Juncker, A. S., Rasmussen, S., Li, J., Sunagawa, S., … MetaHIT Consortium. (2014). Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes. Nature Biotechnology, 32(8), 822–828. https://doi.org/10.1038/nbt.2939

[^2]: Mitchell, A. L., Almeida, A., Beracochea, M., Boland, M., Burgin, J., Cochrane, G., … Finn, R. D. (2019). MGnify: the microbiome analysis resource in 2020. Nucleic Acids Research. https://doi.org/10.1093/nar/gkz1035


Downloading the uhgp-90 database can take a while so expect this demo to be longer than usual.


``` shell
# Go to this repo directory
cd kaamer-demo

# Download prebuilt kAAmer database of Mgnify for gut metagenome profiling
wget https://bacteriapps.genome.ulaval.ca/dbdwl/uhgp-90-2019_09.kaamer.tgz -O data/uhgp-90-2019_09.kaamer.tgz
tar xvf data/uhgp-90-2019_09.kaamer.tgz -C data/

# Download demo gut metagenome fastq 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR209/ERR209293/ERR209293_1.fastq.gz -O data/mg-read.fastq.gz

# Start the Mgnify Human Gut kaamer database...
kaamer-db -server -d data/uhgp-90-2019_09.kaamer &

# Wait for the database to be fully opened ...
kaamer -search -t fastq -i data/mg-reads.fastq.gz -o data/mg-reads.ann.tsv -ann

# Analyse 


```
 
