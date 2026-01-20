# 数据获取
从文章中得知测序所得到的数据的 GEO 编号为 GSE119044  
访问[GEO数据库](https://www.ncbi.nlm.nih.gov/geo/)，得到所需要的 RNA-Seq 数据的 SRA 号为：SRR7753893、SRR7753894、SRR7753895和SRR7753896，均为双端测序  
利用 ENA 数据库，得到所有数据的 FTP 网址，通过 aria2 下载 fastq 文件  

```bash
mkdir -p ~/MC-LR/sequence
cd ~/MC-LR/sequence
aria2c -c -d . -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/003/SRR7753893/SRR7753893_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/003/SRR7753893/SRR7753893_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/004/SRR7753894/SRR7753894_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/004/SRR7753894/SRR7753894_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/005/SRR7753895/SRR7753895_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/005/SRR7753895/SRR7753895_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/006/SRR7753896/SRR7753896_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/006/SRR7753896/SRR7753896_2.fastq.gz
#通过 md5 码查看数据完整性  
md5sum SRR7753893_1.fastq.gz SRR7753893_2.fastq.gz SRR7753894_1.fastq.gz SRR7753894_2.fastq.gz SRR7753895_1.fastq.gz SRR7753895_2.fastq.gz SRR7753896_1.fastq.gz SRR7753896_2.fastq.gz
```