#Variant calling projects, follow steps in Dr Tamer Mansour course NGS2: https://github.com/drtamermansour/nu-ngs02
#Trick on the installation of VEP

#install git
conda install git
#install perl-dbi
conda install -c bioconda perl-dbi

git clone https://github.com/Ensembl/ensembl-vep.git 

cd ensembl-vep
git pull
git checkout release/87
perl INSTALL.pl --NO_HTSLIB

#install perl-mysql
sudo apt-get install libdbd-mysql-perl
OR
conda install -c bioconda perl-dbd-mysql 
 
#get the variant file
perl vep.pl -i BD143_TGACCA_L005.vcf --[cached file for homosapiens]


Sources:

- https://www.biostars.org/p/182780/
-https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html
-https://www.ensembl.org/info/docs/tools/vep/script/index.html
