# 数据获取
从文章中得知测序所得到的数据的 GEO 编号为 GSE119044  
访问[GEO数据库](https://www.ncbi.nlm.nih.gov/geo/)，得到所需要的 RNA-Seq 数据的 SRA 号为：SRR7753893、SRR7753894、SRR7753895和SRR7753896，均为双端测序  
利用[ENA](https://www.ebi.ac.uk/ena/browser/home)数据库，得到所有数据的 FTP 网址，通过 aria2 下载 fastq 文件  

```bash
mkdir -p ~/MC-LR/sequence
cd ~/MC-LR/sequence
aria2c -c -d . -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/003/SRR7753893/SRR7753893_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/003/SRR7753893/SRR7753893_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/004/SRR7753894/SRR7753894_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/004/SRR7753894/SRR7753894_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/005/SRR7753895/SRR7753895_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/005/SRR7753895/SRR7753895_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/006/SRR7753896/SRR7753896_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/006/SRR7753896/SRR7753896_2.fastq.gz
# 通过 md5 码查看数据完整性  
md5sum SRR7753893_1.fastq.gz SRR7753893_2.fastq.gz SRR7753894_1.fastq.gz SRR7753894_2.fastq.gz SRR7753895_1.fastq.gz SRR7753895_2.fastq.gz SRR7753896_1.fastq.gz SRR7753896_2.fastq.gz
```

# 获取参考基因组与注释文件
从[Ensembl 数据库](https://www.ensembl.org/index.html?redirect=no)获取小鼠`Mus musculus`的参考基因组序列  

```bash
mkdir -p ~/MC-LR/genome
cd ~/MC-LR/genome
wget https://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
gzip -d Mus_musculus.GRCm38.dna.toplevel.fa.gz
# 重命名
mv Mus_musculus.GRCm38.dna.toplevel.fa Mus_GRCm38.fa
# 去除描述信息
mv Mus_GRCm38.fa Mus_GRCm38.raw.fa
cat Mus_GRCm38.raw.fa | perl -n -e 'if(m/^>(.+?)(?:\s|$)/){ print ">$1\n";}else{print}' > Mus_GRCm38.fa
rm Mus_GRCm38.raw.fa
```

获取注释信息  

```bash
mkdir -p ~/MC-LR/annotation
cd ~/MC-LR/annotation
wget https://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz
gzip -d Mus_musculus.GRCm38.92.gtf.gz
mv Mus_musculus.GRCm38.92.gtf Mus_GRCm38.gtf
```
# 质量控制与修剪
## 质量评估

```bash
mkdir -p ~/MC-LR/output/fastqc
cd ~/MC-LR/sequence
fastqc -t 4 -o ../output/fastqc *.gz
# 合并
cd ~/MC-LR/output/fastqc
multiqc .
```

## 去除接头和低质量序列
