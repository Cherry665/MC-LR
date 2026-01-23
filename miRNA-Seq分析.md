# 数据下载
从文章中得知测序所得到的数据的 GEO 编号为 GSE119044  
访问[GEO数据库](https://www.ncbi.nlm.nih.gov/geo/)，得到所需要的 RNA-Seq 数据的 SRA 号为：SRR7753897、SRR7753898、SRR7753899和SRR7753900  
利用[ENA](https://www.ebi.ac.uk/ena/browser/home)数据库，得到所有数据的 FTP 网址，通过 aria2 下载 fastq 文件  

```bash
mkdir -p ~/MC-LR/miRNA-Seq/sequence
cd ~/MC-LR/miRNA-Seq/sequence
aria2c -c -d . -Z ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/007/SRR7753897/SRR7753897.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/008/SRR7753898/SRR7753898.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/009/SRR7753899/SRR7753899.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR775/000/SRR7753900/SRR7753900.fastq.gz
# 通过 md5 码查看数据完整性  
md5sum SRR7753897.fastq.gz SRR7753898.fastq.gz SRR7753899.fastq.gz SRR7753900.fastq.gz
```

# miRNA参考序列下载
从[miRBase](https://www.mirbase.org)下载所有物种成熟 miRNA 序列，提取小鼠的 miRNA 序列，利用`bowtie`构建比对索引  

```bash
mkdir -p ~/MC-LR/miRNA-Seq/miRBase
cd ~/MC-LR/miRNA-Seq/miRBase
wget https://www.mirbase.org/download/mature.fa
# 提取小鼠的 miRNA 序列，并完成 U → T 转换
# （在 miRBase 数据库中使用 U 来代表尿嘧啶，以反映miRNA作为RNA分子的原始化学状态，但后续分析软件只使用 A, T, C, G, N 这五个字符）
grep -A 1 '^>mmu-' mature.fa | sed '/^[^>]/ y/Uu/Tt/' > mmu_mature.fa
# 构建比对索引
# （生成4个核心索引文件和6个反向索引文件）
bowtie-build mmu_mature.fa mmu_mature
mkdir ./index
mv *.ebwt ./index
```

# 去除接头和低质量序列
