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
从[miRBase](https://www.mirbase.org)下载所有物种成熟 miRNA 序列(与文章中相同使用了 release 21 版本)，提取小鼠的 miRNA 序列，利用`bowtie`构建比对索引  
mature.fa: 成熟miRNA序列，用于定量；hairpin.fa: 前体序列，用于新miRNA预测  

```bash
mkdir -p ~/MC-LR/miRNA-Seq/miRBase
cd ~/MC-LR/miRNA-Seq/miRBase
wget https://www.mirbase.org/download_version_files/21/mature.fa
# 提取小鼠的 miRNA 序列，并完成 U → T 转换
# （在 miRBase 数据库中使用 U 来代表尿嘧啶，以反映miRNA作为RNA分子的原始化学状态，但后续分析软件只使用 A, T, C, G, N 这五个字符）
grep -A 1 '^>mmu-' mature.fa | sed '/^[^>]/ y/Uu/Tt/' > mmu_mature.fa
# 构建比对索引
# （生成4个核心索引文件和6个反向索引文件）
bowtie-build mmu_mature.fa mmu_mature
rm mature.fa
```

# 质量控制
## 质量评估
```bash
cd ~/MC-LR/miRNA-Seq/sequence
mkdir -p ../output/fastqc
fastqc -t 4 -o ../output/fastqc *.gz
# 合并
cd ~/MC-LR/miRNA-Seq/output/fastqc
multiqc .
```
可以看出测序质量还是可以的，但 3' 端存在接头序列需要去除  
[所有的测序文件的质量](./images/所有测序文件质量.png)  
[平均质量值的read的数量](./images/平均质量值的read的数量.png)  
[接头情况](./images/接头情况.png)  

## 去除接头和低质量序列
miRNA 的一般用`cutadapt`,同时去掉 reads 中的接头，低质量的 reads 以及过长过短的 reads  
--quality-base=33：指定Phred质量分数编码为33  
--minimum-length=18 和 --maximum-length=30：修剪后只保留18-30nt长度的序列  
-q 20：将碱基质量阈值设为20  
--discard-untrimmed：丢弃所有未检测到接头的读段  

```bash
cd ~/MC-LR/miRNA-Seq/sequence
mkdir -p ../output/adapter
for i in $(ls *.fastq.gz);do
    cutadapt -a TGGAATTCTCGGGTGCCAAGG \
    --quality-base=33 \
    --minimum-length=18 \
    --maximum-length=30 \
    -q 20 \
    --discard-untrimmed \
    -o ../output/adapter/${i}  ${i}
done
```

## 再次查看质量情况
```bash
cd ~/MC-LR/miRNA-Seq/output/adapter
mkdir -p ../fastqc_adapter
fastqc -t 4 -o ../fastqc_adapter *.gz
cd ../fastqc_adapter
multiqc .
```

# 序列比对
-n：允许错配的数量  
-m：允许比对到参考序列的最多条数  
--best --strata：生成的sam文件只显示最佳的 map 结果  
-S：输出结果文件为 sam 格式  
```bash
mkdir -p ~/MC-LR/miRNA-Seq/output/align
cd ~/MC-LR/miRNA-Seq/output/adapter
parallel -k -j 4 "
    bowtie -n 2 -m 10 --best --strata  -x ../../miRBase/mmu_mature {1}.fastq.gz  -S ../align/{1}.sam 2>../align/{1}.log
" ::: $(ls *.fastq.gz | perl -p -e 's/\.fastq\.gz$//')
# 转换为 BAM 格式
cd ~/MC-LR/miRNA-Seq/output/align
parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
# 统计 BAM 文件中比对到各个参考序列的 reads 数量  
# 生成的.txt 文件每列的含义分别为：miRNA_ID  长度  比对到该序列的reads数  未比对到该序列的reads数
parallel -k -j 4 "
    samtools idxstats {1}.sort.bam > {1}.txt
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
# 合并表达矩阵
paste *.txt |cut -f 1,3,7,11,15 > mmu.txt
# 加上列名
echo -e "miRNA\tSRR7753897\tSRR7753898\tSRR7753899\tSRR7753900" | cat - mmu.txt > mmu.mature.txt
```  

# 差异表达分析
无生物学重复使用 DEGseq，有生物学重复使用 DESeq2  

```R

```
