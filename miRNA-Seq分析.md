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
wget https://www.mirbase.org/download_version_files/21/hairpin.fa
# 提取小鼠的 miRNA 序列，并完成 U → T 转换
# （在 miRBase 数据库中使用 U 来代表尿嘧啶，以反映miRNA作为RNA分子的原始化学状态，但后续分析软件只使用 A, T, C, G, N 这五个字符）
# mature.fa
grep '^>mmu-' mature.fa | sed 's/^>//' > mmu_mature_names.txt
faops some mature.fa mmu_mature_names.txt mmu_mature.fa
sed -i '/^[^>]/ y/Uu/Tt/' mmu_mature.fa
# hairpin.fa
grep '^>mmu-' hairpin.fa | sed 's/^>//' > mmu_hairpin_names.txt
faops some hairpin.fa mmu_hairpin_names.txt mmu_hairpin.fa
sed -i '/^[^>]/ y/Uu/Tt/' mmu_hairpin.fa
# 构建比对索引
# （生成4个核心索引文件和6个反向索引文件）
bowtie-build mmu_mature.fa mmu_mature
bowtie-build mmu_hairpin.fa mmu_hairpin
rm mature.fa hairpin.fa
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
## mature miRNA
-n：允许错配的数量  
-m：允许比对到参考序列的最多条数  
--best --strata：生成的sam文件只显示最佳的 map 结果  
-S：输出结果文件为 sam 格式  

```bash
mkdir -p ~/MC-LR/miRNA-Seq/output/align/mature
cd ~/MC-LR/miRNA-Seq/output/adapter
parallel -k -j 4 "
    bowtie -n 2 -m 10 --best --strata  -x ../../miRBase/mmu_mature {1}.fastq.gz  -S ../align/mature/{1}.sam 2>../align/mature/{1}.log
" ::: $(ls *.fastq.gz | perl -p -e 's/\.fastq\.gz$//')
# 转换为 BAM 格式
cd ~/MC-LR/miRNA-Seq/output/align/mature
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

## hairpin miRNA
```bash
mkdir -p ~/MC-LR/miRNA-Seq/output/align/hairpin
cd ~/MC-LR/miRNA-Seq/output/adapter
parallel -k -j 4 "
    bowtie -n 2 -m 10 --best --strata  -x ../../miRBase/mmu_hairpin {1}.fastq.gz  -S ../align/hairpin/{1}.sam 2>../align/hairpin/{1}.log
" ::: $(ls *.fastq.gz | perl -p -e 's/\.fastq\.gz$//')
# 转换为 BAM 格式
cd ~/MC-LR/miRNA-Seq/output/align/hairpin
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
echo -e "miRNA\tSRR7753897\tSRR7753898\tSRR7753899\tSRR7753900" | cat - mmu.txt > mmu.hairpin.txt
```

# 差异表达分析
无生物学重复使用`DEGseq`，有生物学重复使用`DESeq2`(以下使用的是`DESeq2`)  

```R
# 读取.txt 文件并创建数据框，将第一行作为列名，将第一列作为行名
dataframe <- read.table("mmu.mature.txt", header=TRUE, row.names = 1)
# 去除低表达基因（基因在所有样本中的表达量总和不为0）
countdata <- dataframe[rowSums(dataframe) > 0,]
# 安装加载所需要的 R 包
# 加载包
library(DESeq2)
library(biomaRt)
# 构建对象
sample_names <- c("SRR7753897", "SRR7753898", "SRR7753899", "SRR7753900")
condition <- c("control", "control", "MC_LR", "MC_LR")
coldata <- data.frame(row.names = sample_names, condition = condition)
# 调整数据顺序
countdata <- countdata[row.names(coldata)]
# 构建 dds （design 后面的词需和 coldata 的列名相同）
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ condition)
# 设置因子水平的顺序，让 DESeq2 在比较时以 control 组作为对照组
dds$condition <- factor(dds$condition, levels = c("control", "MC_LR"))
# 对差异基因进行分析，并进行统计检验
dds <- DESeq(dds)
res <- results(dds)
# 提取差异基因（使用论文阈值）
# sig_genes <- subset(res, abs(log2FoldChange) > 1 & pvalue < 0.05)
# 绘制火山图（利用ggplot2）
library(ggplot2)
log2FC_threshold <- 1
pvalue_threshold <- 0.05
log10p_threshold <- -log10(pvalue_threshold)
plot_data <- as.data.frame(res)
plot_data$group <- "NS" 
plot_data$group[which(plot_data$pvalue < 0.05 & plot_data$log2FoldChange > 1)] <- "Up"
plot_data$group[which(plot_data$pvalue < 0.05 & plot_data$log2FoldChange < -1)] <- "Down"
p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), 
             linetype = "dashed", color = "grey", alpha = 0.5, linewidth = 0.7) +
  geom_hline(yintercept = log10p_threshold, 
             linetype = "dashed", color = "grey", alpha = 0.5, linewidth = 0.7) +
  geom_point(aes(color = group), alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal()
# 保存图片为 pdf 文件
ggsave("miRNA差异表达分析.pdf", plot = p, width = 8, height = 6, dpi = 300)
```