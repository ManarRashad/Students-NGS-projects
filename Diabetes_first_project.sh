## important steps before download data, Dr Tamer Mansour: 
#https://github.com/drtamermansour/nu-ngs01
#install Bioconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  ## follow the prompt. Keep pressing ENTER or responding by yes when needed :)
# restart the terminal
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
#Download sra-tools:
conda install -c bioconda sra-tools
#Download entrez-direct
conda install -c bioconda entrez-direct

##Download data
esearch -db sra -query [ID_of _project] | efetch -format runinfo | cut -d "," -f 1 | grep SRR | xargs fastq-dump -X 10000 --skip-technical --read-filter pass –dumpbase --gzip

##install Fastqc and apply it
#install it
conda install -c bioconda fastqc 
conda install -c bioconda multiqc
#apply it
for f in *.fq.gz;do fastqc -t 1 -f fastq -noextract $f;done

## important steps before apply kallisto

# Download human transcriptome
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz

gunzip gencode.v29.pc_transcripts.fa.gz

# Download the Transcriptome Annotation File
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

gunzip gencode.v29.annotation.gtf.gz

# Select the transcripts of Chr22

READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS
    do  
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    done

## install Kallisto and apply it
#install it
conda install -c bioconda -y kallisto

# apply it
#indexing
kallisto index -i human_pc.idx -k 25 gencode.v29.pc_transcripts.chr22.simplified.fa

#pseudoalignment
for value in *fastq; do kallisto quant -i human_pc.idx -o $value.human_pc_bam $value --pseudobam --single -l 200 -s 20; done

#convert the kallisto output to Deseq input and run the deseq2 then perform functional enrichment analysis using Gprofiler
Rscript taximport1.r
