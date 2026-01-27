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

## 获取注释信息  

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

```bash
mkdir -p ~/MC-LR/output/trim
cd ~/MC-LR/sequence
parallel -j 4 "
  java -jar ~/biosoft/Trimmomatic-0.38/trimmomatic-0.38.jar \
    PE -phred33 {1} {= s/_1\.fastq\.gz/_2.fastq.gz/ =} \
    ~/MC-LR/output/trim/{= s/_1\.fastq\.gz/_1_paired.fq.gz/ =} \
    ~/MC-LR/output/trim/{= s/_1\.fastq\.gz/_1_unpaired.fq.gz/ =} \
    ~/MC-LR/output/trim/{= s/_1\.fastq\.gz/_2_paired.fq.gz/ =} \
    ~/MC-LR/output/trim/{= s/_1\.fastq\.gz/_2_unpaired.fq.gz/ =} \
    ILLUMINACLIP:$HOME/biosoft/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *_1.fastq.gz)
rm ../output/trim/*_unpaired.fq.gz
```

# 再次查看质量情况
```bash
mkdir -p ~/MC-LR/output/fastqc_trim
cd ~/MC-LR/output/trim
fastqc -t 4 -o ../fastqc_trim *.gz
cd ../fastqc_trim
multiqc .
```

# 序列比对
## 下载索引
```bash
mkdir -p ~/MC-LR/genome/index
cd ~/MC-LR/genome/index
wget https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download
# 下载下来的 download 文件是.tar.gz 文件  
tar -xzvf download
```

## 序列比对
```bash
mkdir -p ~/MC-LR/output/align
cd ~/MC-LR/output/trim
parallel -k -j 4 "
    hisat2 -t -x ../../genome/index/grcm38/genome \
      -1 {1}_1_paired.fq.gz -2 {1}_2_paired.fq.gz \
      -S ../align/{1}.sam \
      --summary-file ../align/{1}.summary \
      2>../align/{1}.log
" ::: $(ls *_1_paired.fq.gz | perl -p -e 's/_1_paired\.fq\.gz$//')

# 将.sam 文件压缩为.bam 文件,输出排序后的.bam 文件和对应的索引文件  
cd ~/MC-LR/output/align
parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

# 表达量统计
```bash
mkdir -p ~/MC-LR/output/htseq
cd ~/MC-LR/output/align
# 将.bam 文件改为按 name 进行排序  
parallel -j 4 "
    samtools sort -n -@ 4 -o {1}.name_sorted.bam {1}.sort.bam
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')
# 进行表达量统计  
parallel -j 4 "
    htseq-count -s no -f bam -r name {1}.name_sorted.bam ../../annotation/Mus_GRCm38.gtf \
      >../htseq/{1}.count  2>../htseq/{1}.log
" ::: $(ls *.name_sorted.bam | perl -p -e 's/\.name_sorted\.bam$//')
```

# 合并表达矩阵
```R
rm(list=ls())
setwd("~/MC-LR/output/htseq")
# 得到文件样本编号
files <- list.files(".", "*.count")
f_lists <- list()
for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[prefix]] = i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  data[[count]] <- a
}
# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]],by="gene_id")
}
write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```

# 差异表达分析
## 分析前处理
```R
# 读取.csv 文件并创建数据框，将第一行作为列名，将第一列作为行名  
dataframe <- read.csv("merge.csv", header=TRUE, row.names = 1)
# 去除前5行描述行（不包括列名）
countdata <- dataframe[-(1:5),]
# 去除 ID （此处不用，基因 ID 带版本号时需要）
row_names <- row.names(countdata)
name_replace <- gsub("\\.\\w+","", row.names(countdata))
row.names(countdata) <- name_replace
# 去除低表达基因（基因在所有样本中的表达量总和不为0）
countdata <- countdata[rowSums(countdata) > 0,]
```

## 差异分析
```R
# 安装加载所需要的 R 包
# 下载小鼠基因的注释信息数据库 org.Mm.eg.db
BiocManager::install("org.Mm.eg.db")
# 加载包
library(DESeq2)
library(biomaRt)
library(org.Mm.eg.db)
# 构建对象
sample_names <- c("SRR7753893", "SRR7753894", "SRR7753895", "SRR7753896")
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
plot_data <- as.data.frame(res)
plot_data$group <- "NS" 
plot_data$group[which(plot_data$pvalue < 0.05 & plot_data$log2FoldChange > 1)] <- "Up"
plot_data$group[which(plot_data$pvalue < 0.05 & plot_data$log2FoldChange < -1)] <- "Down"
p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = group), alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal()
# 保存图片为 pdf 文件
ggsave("MC-LR差异表达分析.pdf", plot = p, width = 8, height = 6, dpi = 300)
```

# GO 富集分析（利用 ClusterProfiler）
```R
library(clusterProfiler)
# 获取显著上调/下调基因
up_genes <- rownames(res[res$log2FoldChange > 1 & res$pvalue < 0.05, ])

down_genes <- rownames(res[res$log2FoldChange < -1 & res$pvalue < 0.05, ])
# GO富集分析（关注“生物过程BP”）
go_up <- enrichGO(gene = up_genes,
                  OrgDb = org.Mm.eg.db,
                  keyType = 'ENSEMBL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_down <- enrichGO(gene = down_genes,
                  OrgDb = org.Mm.eg.db,
                  keyType = 'ENSEMBL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
# 绘制条形图
a <- barplot(go_up, showCategory=8, title="GO Enrichment (Up-regulated)")
b <- barplot(go_down, showCategory=8, title="GO Enrichment (Down-regulated)")
# 保存
ggsave("MC-LR上调基因 GO 富集分析.pdf", plot = a, width = 8, height = 6, dpi = 300)
ggsave("MC-LR下调基因 GO 富集分析.pdf", plot = b, width = 8, height = 6, dpi = 300)

# 保存数据
write.table(res, "./DE_result.tsv", sep="\t", quote = FALSE)
write.table(up_genes, "./up_genes.tsv", sep="\t", quote = FALSE)
write.table(down_genes, "./down_genes.tsv", sep="\t", quote = FALSE)
```