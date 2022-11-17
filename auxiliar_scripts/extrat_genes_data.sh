curl -O http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.104.chr.gtf.gz
sed  's/;/\t/g' Homo_sapiens.GRCh38.104.chr.gtf | less > genes_hg38.txt