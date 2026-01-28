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
# mature.fa
grep '^>mmu-' mature.fa | sed 's/^>//' > mmu_mature_names.txt
faops some mature.fa mmu_mature_names.txt mmu_mature.fa
sed -i '/^[^>]/ y/Uu/Tt/' mmu_mature.fa
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
## mature miRNA
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
## DESeq2
```R
# 读取.txt 文件并创建数据框，将第一行作为列名，将第一列作为行名
dataframe <- read.table("mmu.mature.txt", header=TRUE, row.names = 1)
# 去除低表达序列（基因在所有样本中的表达量总和不为0）
countdata <- dataframe[rowSums(dataframe) > 0,]
# 安装加载所需要的 R 包
# 加载包
library(DESeq2)
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
# 对差异表达序列进行分析，并进行统计检验
dds <- DESeq(dds)
res <- results(dds)
# 提取差异表达序列（使用论文阈值）
up <- rownames(res[res$log2FoldChange > 1 & res$pvalue < 0.05, ])
down <- rownames(res[res$log2FoldChange < -1 & res$pvalue < 0.05, ])
# 保存数据
write.table(res, "./miRNA_DE_result.tsv", sep="\t", quote = FALSE)
write.table(up, "./miRNA_DE_up.tsv", sep="\t", quote = FALSE)
write.table(down, "./miRNA_DE_down.tsv", sep="\t", quote = FALSE)
```
## edgeR
```R
library(edgeR)

```

# 可视化（利用ggplot2）
```R
library(ggplot2)
# 绘制火山图  
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
ggsave("miRNA差异表达分析DESeq2.pdf", plot = p, width = 8, height = 6, dpi = 300)
```

# 靶基因预测
## 下载miRanda
```bash
conda install bioconda::miranda
```

## 使用 miRanda 预测靶基因
```bash
# 激活conda环境  
conda activate  

# 获取上调/下调序列的.fa 文件  
mkdir -p ~/MC-LR/miRNA-Seq/output/miRanda
cd ~/MC-LR/miRNA-Seq/miRBase
faops some mmu_mature.fa ../output/align/miRNA_DE_up.tsv ../output/miRanda/mmu_up.fa
faops some mmu_mature.fa ../output/align/miRNA_DE_down.tsv ../output/miRanda/mmu_down.fa
```
3'UTR 序列文件是一个专门包含基因3'非翻译区 DNA 序列的文本文件。绝大多数 miRNA 通过与靶基因 mRNA 的 3'UTR 区域碱基互补配对，从而抑制翻译或导致 mRNA 降解  

```R
# 下载小鼠的 3’UTR 序列文件
cd ~/MC-LR/miRNA-Seq/output/miRanda
library(biomaRt)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
## 建立与 Ensembl 数据库的连接，指定使用 Ensembl 的“Mart”服务，并选择小鼠基因数据集
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
## 获得小鼠基因转录本ID
ensembl_ID <- toTable(org.Mm.egENSEMBLTRANS)
## 获得 3’UTR 序列文件
utr <- getSequence(id=ensembl_ID$trans_id, type="ensembl_transcript_id", seqType='3utr', mart=ensembl)
## 获得 3’UTR 序列的.fa文件
## utr 中第1列为序列，第2列为转录本ID
outfile <- file("mmu-3utr.fa", "w")
for (i in 1:nrow(utr)) {
  if(grepl("Sequence unavailable", utr[i, 1])) {
    next
  }
  h = paste(c(">", utr[i,2]), collapse="")
  writeLines(h, outfile)
  writeLines(utr[i,1], outfile)
}
close(outfile)
```
```bash
# 靶基因预测
# -sc：指定序列比对打分的阈值，小于该阈值的结合位点会被过滤，默认为140
## 获取表达上调 miRNA 的靶基因转录本ID
miranda mmu_up.fa mmu-3utr.fa -out up_target.txt
grep '^>' up_target.txt | grep -v '>>' > up_result.txt
cut -f 2 up_result.txt > up_genes.txt
## 获取表达下调 miRNA 的靶基因转录本ID
miranda mmu_down.fa mmu-3utr.fa -out down_target.txt
grep '^>' down_target.txt | grep -v '>>' > down_result.txt
cut -f 2 down_result.txt > down_genes.txt
```

# GO 富集分析
```R
library(org.Mm.eg.db)
library(clusterProfiler)
# 将 Ensembl 转录本ID转换为 ENTREZID
# 可能会出现同一个基因产生多个不同的转录本或一个转录本关联到多个可能的基因上的情况（比例少不做去除）
up_targetname= read.csv("up_genes.txt",head=F)
up_gene.df <- bitr(up_targetname$V1, fromType = "ENSEMBLTRANS",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

down_targetname= read.csv("down_genes.txt",head=F)
down_gene.df <- bitr(down_targetname$V1, fromType = "ENSEMBLTRANS",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
# GO 富集分析（关注“生物过程BP”）
up_ego <- enrichGO(gene = up_gene.df$ENTREZID,
                OrgDb = org.Mm.eg.db,
                keyType = 'ENTREZID',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

down_ego <- enrichGO(gene = down_gene.df$ENTREZID,
                OrgDb = org.Mm.eg.db,
                keyType = 'ENTREZID',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
# 可视化
library(ggplot2)
a <- barplot(up_ego, showCategory=8, title="GO Enrichment (Up-regulated)")
B <- barplot(down_ego, showCategory=8, title="GO Enrichment (Down-regulated)")

```
