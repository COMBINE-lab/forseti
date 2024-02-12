# Forseti
### This repo is for the Forseti Python package. For the code and scripts used for reproducing the results in the Forseti paper, please refer to [forseti-experiments](https://github.com/COMBINE-lab/forseti-experiments).

**_Forseti_** is a predictive model to probabilistically assign a splicing status to scRNA-seq reads. Our model has two key components: first, we train a binding affinity model to assign a probability that a given transcriptomic site is used in fragment generation; second, we fit a robust fragment length distribution model that generalizes well across datasets deriving from different species and tissue types.  

**_Forseti_** combines these two trained models to predict the splicing status of the molecule of origin of reads by scoring putative fragments that associate each alignment of sequenced reads with proximate potential priming sites. 

In this tutorial, we will show you how to predict the splicing status of scRNA-seq for a standard scRNA-Seq dataset using _Forseti_. 

## Set up and download _Forseti_
Before start, you can run the following commands to create a working directory and download the _Forseti_ repo from Github. We alias this directory and use the alias below so that you can easily set it to something else and still copy and paste the later shell commands.  

```bash
# Run in shell
# Set up the working directory
export $FST_SAMPLE_DIR=$PWD/forseti_sample_dir
mkdir -p $FST_SAMPLE_DIR

# Download the Forseti package
cd $FST_SAMPLE_DIR
git clone https://github.com/COMBINE-lab/forseti.git .
```
## Download sample data and references 
For this tutorial, we use sample dataset [_PBMC1k (v3) healthy donor samples_](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0) from 10X website.
```bash
# Download the FASTQ files 
FASTQ_DIR="$FST_SAMPLE_DIR/data/pbmc_1k_v3_fastqs"
mkdir -p $FASTQ_DIR

wget -qO- https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar | tar xf - --strip-components=1 -C $FASTQ_DIR

# Download human reference genome build and gene annotations
export REF_DIR="$FST_SAMPLE_DIR/data/refdata-gex-GRCh38-2020-A"
mkdir -p $REF_DIR

wget -qO- https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz | tar xzf - --strip-components=1 -C $REF_DIR
```

## Set up environments
> Notes: To handle the potential conflicts for environments having both R and python, we have two seperate envs:
* _build_spliceu_r_ for building spliceu index.
* _forseti_test_ for the rest of the tutorial.
```bash
# Create the env for building spliceu index
conda env create -f $FST_SAMPLE_DIR/forseti/envs/build_spliceu_r.yml
# Activate the env
conda activate build_spliceu_r
```
## 1- Build the _spliceu_ index
The R script `build_spliceu.R` will build a _spliceu_ (spliced + unspliced), an augmented transcriptome reference consisting of a augmented gene annotation GTF file, a transcriptome FASTA file, and a transcript-to-gene mapping TSV file. These files will be used in the later steps to align reads and predict the splicing status. 

**Input**: A standard reference, including a gene annotation file(.gtf) and the corresponding genome build (.fasta).  
**Output**: The _spliceu_ augmented reference: spliceu.gtf, spliceu.fa, t2g_3col.tsv   

```bash
# Run in shell
export IDX_DIR="$FST_SAMPLE_DIR/data/spliceu_refs"
Rscript $FST_SAMPLE_DIR/forseti/preprocess/build_spliceu.R $REF_DIR/fasta/genome.fa $REF_DIR/genes/genes.gtf $IDX_DIR 12
mv $REF_DIR/fasta/genome.fa $IDX_DIR/genome.fa
mv $REF_DIR/genes/genes.gtf $IDX_DIR/genes.gtf
```

```bash
# When you completed the step-1, deactivate its env with this command
conda deactivate 
```
### Then, you can use the yml file to create an environment for _Forseti_. 

```bash
conda env create -f $FST_SAMPLE_DIR/forseti/envs/forseti.yml
conda activate forseti_test
```

## 2- Align with STAR 
With the following commmands, generate genome indices for STAR.
```bash
# Run in shell
STAR --runThreadN 12 \
--runMode genomeGenerate --genomeDir $IDX_DIR \
--genomeFastaFiles $IDX_DIR/genome.fa \
--sjdbGTFfile $IDX_DIR/spliceu.gtf \
--sjdbOverhang 100
```
Run the scripts below to align the downloaded PBMC1k sample data with STAR. 

**Input**: reads in FATSQ files, spliceu indices  
**Output**: uniquely mapped reads (.bam)
```bash
# Align the reads to reference
STAR_IDX_REFS="$FST_SAMPLE_DIR/data/star_spliceu_refs"

export STAR_out=$FST_SAMPLE_DIR/STAR_out
mkdir -p $STAR_out

STAR --runThreadN 12 \
    --readFilesIn $FASTQ_DIR/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,$FASTQ_DIR/pbmc_1k_v3_S1_L002_R2_001.fastq.gz $FASTQ_DIR/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,$FASTQ_DIR/pbmc_1k_v3_S1_L002_R1_001.fastq.gz\
    --readFilesCommand zcat \
    --genomeDir $STAR_IDX_DIR \
    --outFileNamePrefix $STAR_out \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --soloType CB_UMI_Simple \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist $FST_SAMPLE_DIR/data/whitelist/3M-february-2018.txt \
    --soloFeatures GeneFull \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --limitIObufferSize 50000000 50000000 
```

- `--readFilesIn` passes the path to your fastq file.
- Use `--readFilesCommand zcat` if you have compressed fastq files (i.e. *.gz)
- With `--genomeDir` option, you pass the path to the folder storing the spliceu index. 
- `--outFileNamePrefix` passes the path to your output
- `--outSAMtype BAM SortedByCoordinate` will output the sorted BAM file 
- With `--quantMode TranscriptomeSAM`, STAR will output alignments translated into transcript coordinates, and `--quantTranscriptomeBan Singleend` allows insertions, deletions ans soft-clips in the transcriptomic alignments.
- The PBMC1k data in this tutorial is 10X chemistry v3, so we use `--soloType CB_UMI_Simple` , `--soloUMIlen 12`, `--soloBarcodeReadLength 0` and `--soloCBwhitelist 3M-february-2018.txt`.
- `--soloFeatures GeneFull` outputs exon and intron UMI counts.


> Note: The following section is an OPTINAL step for standard process. **In short**, we extract the top cells to reduce the time for the tutorial. Therefore,
- if you are **playing with the 10X 1kPBMC sample data** and want to save time and check how the predcition model works, you could follow this step. 
- If you are **using a small dataset**, you can skip this section and continue to the next step: _filter the reads have unique map to the genome_.

The 10X 1kPBMC sample data have 200 million algnments when mapping to Transcriptom. To reduce the time taken for completing the tutorial, **we extract the top 200 cells**. The filtering is done based on the gene count matrix provided on 10X website. The top cells defined by cells with at least 3000 genes in count matrix.
If you are playing with the 10X 1kPBMC sample data, you can check the scripts we used to filter the reads, which is saved in `forseti/preprocess/get_top_cell.sh` and `extract_top_cells.py`. The filtered bam file contains 65,784,347 alignments, and will take 3.3 hour in 32threads for the prediciton.


Then, we use samtools to filter the reads have unique map to the genome. 
```bash
# Run in shell
STAR_OUT_DIR="$FST_SAMPLE_DIR/STAR_out"
# filter genome unique mapped reads 
# `-d NH:1` with this flag, we get the uniquely mapped reads(i.e. number of hit is 1).
samtools view -d NH:1 $STAR_OUT_DIR/Aligned.sortedByCoord.out.bam -o $STAR_OUT_DIR/filtered_Aligned.sortedByCoord.out.bam

# extract the list of unique map reads
samtools view $STAR_OUT_DIR/filtered_Aligned.sortedByCoord.out.bam | cut -f1 > $STAR_OUT_DIR/unimap_read_names.txt
# extract the toTx record based on the unique maped read list
samtools view \
    -N $STAR_OUT_DIR/unimap_read_names.txt \
    -o $STAR_OUT_DIR/filtered_Aligned.toTranscriptome.out.bam \
    $STAR_OUT_DIR/Aligned.toTranscriptome.out.bam

# find the reads 100% mapped in a exonic region (i.e. ambiguous reads)
bedtools intersect \
-a $STAR_OUT_DIR/filtered_Aligned.sortedByCoord.out.bam \
-b $IDX_DIR/exon_by_tx.bed \
-f 1.0 \
-wo \
-bed \
> $STAR_OUT_DIR/exonic_reads.bed

# extract the list of exonic reads
cut -f 4 $STAR_OUT_DIR/exonic_reads.bed > $STAR_OUT_DIR/exonic_read_name.txt
# filter the exonic reads in toTranscriptome.bam file
samtools view -N $STAR_OUT_DIR/exonic_read_name.txt -o $STAR_OUT_DIR/exonic_Aligned.toTranscriptome.out.bam $STAR_OUT_DIR/filtered_Aligned.toTranscriptome.out.bam
```
## 3- Apply _Forseti_ prediction model
There are several things you want to know before playing with the model.
> 1. What kind of reads we will cover in this model? 

The model targets reads that map uniquely to a single location on the genome, filtering out any reads that are mapped to multiple genomic positions. Despite their unique mapping, these reads may still present compatibility with transcripts in both splicing statuses—spliced and unspliced—introducing **splicing status ambiguity**. In addition, reads mapped to regions where gene overlap occurs can lead to gene origin ambiguity,leading to the **gene origin ambiguity**. 

> 2. How we classify the reads before we run the predcition model?

Before running our prediction model, we classify uniquely mapped reads into three categories based on their alignments to the transcriptomic reference. These categories are 'U' for unspliced, 'S' for spliced, and 'A' for ambiguous reads. Considering that reads may be compatible with multiple transcripts, we identify these as potential transcripts (abbreviated as 'Tx') of origin for the read. The read is considered unspliced ('U') if all of its potential Tx are unspliced versions; it is classified as spliced ('S') if all potential Tx are spliced versions. 

 Ambiguous reads ('A') are those compatible with both spliced and unspliced transcripts. Our prediction model processes these ambiguous reads and accurately predict splicing status for them. Note, some reads can remain ambiguous after the prediction. 

----------

#### Now, you are ready to use the Forseti model! _Forseti_ is implemented in Python and you could easily import the `forseti` package, then call the master function named `forseti_predictor` for the splicing status prediction.

Overall, you need 6 inputs: 
* 4 files from previous steps
* 2 built-in models 
* 1 parameter setting section. 

Forseti will return 2 objects:
* `forseti_dict` for prediction results
* `predicitons_log` for log information.   

`forseti_dict` is a dictionary of _{ readname : read_predictions }_. The _read_predictions_ was saved with `ForsetiReads` object and `ForsetiPredcitions` object.

### Here is a detailed example for how to use the prediction model:  

### 3.0- Import required packages
> Note: All scripts below (3.0-3.3) will be executed in _Python_. Make sure you created this python file on the top level of the `forseti_sample_dir` folder.
```Python
import os
import pickle
from forseti.prediction_model import forseti_predictor, get_initial_splicing_status
```
### 3.1- Prepare input files.
Specify the path to _spliceu fasta_ file and _t2g_ file you've got from `Step1-Build spliceu index`, and the toTxome bam file (use Trasncriptom as reference) from `Step2-Align with STAR`. Also, the _unimap read names_ derived from the toGenome bam file(use Genome as reference).

```Python
# Pass the paths to input files
forseti_sample_dir= os.getcwd()
idx_dir = os.path.join(forseti_sample_dir, 'data','spliceu_refs')

spliceu_fasta_file = os.path.join(idx_dir, "spliceu.fa")
t2g_file = os.path.join(idx_dir, "t2g_3col.tsv")

reads_txome_bam_file = os.path.join(
    forseti_sample_dir, 'STAR_out', 'Aligned.toTranscriptome.out.bam')
unimap_read_names_file = os.path.join(
    forseti_sample_dir, 'STAR_out','unimap_read_names.txt')
```
We provide 2 pre-built models: the `spline model` for fragment length distribution and the `mlp model` for binding affinity prediction. 
```Python
# Pass the path for spline model and mlp model
spline_model_file = os.path.join(forseti_sample_dir, "forseti","pre_built_models", "spline_model.pkl")
mlp_model_file = os.path.join(forseti_sample_dir, "forseti","pre_built_models","mlp.pkl")
```
Set up parameters if you don't want to use the default values listed below.
```Python
# Parameter setting (default)
snr_min_size = 6 
discount_perc = 1 
polya_tail_len = 200
max_frag_len = 1000
num_threads = 12
```
- With `snr_min_size`, you can define the minimum length of valid adenine-single nucleotide repeat(aSNR). The mlp model will only look for the priming window containing valid aSNR. The default length is 6, suggesting that priming window should contain an aSNR of length at least 6(i.e. 'AAAAAA'). 
- `discount_perc` enables to apply a discount ratio on internal polyA against terminal polyA. The default is 1, which means the intrernal polyA and terminal polyA have same weight and there is no discount.
- `polya_tail_len` is the length of polyA tailed we extend at the end of the transcript. 
- `max_frag_len` represent the maximum range the model would search for the priming window. The default is 1000 bp, we will search the downstream 1000 for sense alignments and upstream 1000 for reverse alignments.
- For this model, `num_threads` is not the bottleneck, just for parallel processing. You can use as many as you have. 


### 3.2- Apply the prediciton function
Yippee 🎉! Now, you have done with preparing all input files, models and the parameter setting. You can directly call the prediction function with the following scripts.
```Python
forseti_dict, predicitons_log = forseti_predictor(reads_txome_bam_file, unimap_read_names_file, spliceu_fasta_file, t2g_file, spline_model_file,mlp_model_file, polya_tail_len, discount_perc, snr_min_size, num_threads)
```
### 3.3- Access to the prediction results

Here, we use different examples to show you how to access to the prediction record for each read.  

First, we get the prediction record from _forseti dict_ by its read name. For each `forseti_read`, you can obtain its read name and boolean for whether it is multi-mapped. 
* In our model, if'_is_multi_mapped_' is TRUE, which means the we got one best gene for the read's origin; while '_is_multi_mapped_' is FALSE, indicates that there are multiple best genes and they got a tie, and all of them will be listed in the _predicitons_.

For the `predictions` :
* _splicing_status_ : 
    - 'U' for 'Unspliced read'
    - 'S' for 'Spliced read'
    - 'A' for 'Ambiguous read'
* _gene_id_ : the gene id of origin
* _orientation_ :
    - '+' for sense 
    - '-' for antisense
* _max_prob_ : represents the maximum probobability when comparing antisense_prob and sense_prob for this gene
* _unsplice_prob_ : probobability of reads is unsplice
* _splice_prob_ : probobability of reads is splice

For reads have unique assigned gene(i.e. forseti_read.is_multi_mapped == True), you can get the only predciiton record by calling '.get_get_predictions()', and then access the specific prediction record values. Below, we show what does the prediciton results look like in 4 examples. Reads are assigned to S(Spliced), U(Unspliced), A(Ambiguous), and multi-mapped reads.
```Python
# example 1
rname = 'A00228:279:HFWFVDMXX:1:1109:30454:13463'
forseti_read = forseti_dict[rname]
print(forseti_read.read_name) # A00228:279:HFWFVDMXX:1:1109:30454:13463
print(forseti_read.is_multi_mapped) # False
prediction = forseti_read.get_predictions()
print(prediction.splicing_status) # S ,indicates spliced
print(prediction.gene_id) # ENSG00000115282
print(prediction.orientation) # +, indicates sense strand
print(prediction.max_prob) # 6.176078485466973e-05
print(prediction.unsplice_prob) # 0.0
print(prediction.splice_prob) # 1.0

# example 2
rname = 'A00228:279:HFWFVDMXX:1:1117:25518:32972'
forseti_read = forseti_dict[rname]
print(forseti_read.read_name) # 'A00228:279:HFWFVDMXX:1:1117:25518:32972'
print(forseti_read.is_multi_mapped) # False
prediction = forseti_read.get_predictions()
print(prediction.splicing_status) # U ,indicates unspliced
print(prediction.gene_id) # ENSG00000117616
print(prediction.orientation) # + , indicates sense strand
print(prediction.max_prob) # 0.006903024469278407
print(prediction.unsplice_prob) # 0.6380657291449385
print(prediction.splice_prob) # 0.36193427085506147

# example 3 for Ambiguous read
rname = 'A00228:279:HFWFVDMXX:1:1101:6994:26256'
forseti_read = forseti_dict[rname]
print(forseti_read.read_name) # 'A00228:279:HFWFVDMXX:1:1117:25518:32972'
print(forseti_read.is_multi_mapped) # False
prediction = forseti_read.get_predictions()
print(prediction.splicing_status) # A ,indicates ambiguous
print(prediction.gene_id) # ENSG00000026025
print(prediction.orientation) # + , indicates sense strand
print(prediction.max_prob) # 0.0003280696383465825
print(prediction.unsplice_prob) # 0.5
print(prediction.splice_prob) # 0.5

# example 4 for how to access the predcition for multimapped reads
rname = 'A00228:279:HFWFVDMXX:1:1101:25373:1783'
forseti_read = forseti_dict[rname]
print(forseti_read.read_name) 
print(forseti_read.is_multi_mapped) 

prediction_lst = forseti_read.get_predictions()
for prediction in prediction_lst:
    print(prediction.gene_id) 
    print(prediction.orientation)
    print(prediction.splicing_status)
    print(prediction.max_prob)
    print(prediction.unsplice_prob)
    print(prediction.splice_prob)

'''
For the multi-mapped reads. it will return the list of predcitions:

A00228:279:HFWFVDMXX:1:1109:17400:19820
True
ENSG00000173821
+
A
0.0068454212748862405
0.5
0.5
ENSG00000263069
-
U
0.0068454212748862405
1.0
0.0
'''

```
Then, the log shows the summarized information of the prediction:
* _read_count_ : the total number of reads.
* _final_multi_gene_count_ : In the final prediction, we count how many reads are multi-gene mapped.
* _rescued_multigene_count_ : This count includes reads that initially align with multiple genes, but have a determinative gene origin by comparing the prediction scores.
* _final_ambiguous_count_ : For reads that are ultimately mapped to a unique gene (originally mapped to a unique gene or rescued), if their final splicing status remains ambiguous, they are included in the 'final_ambiguous_count'.
```Python
# Check the summarized log info
print(f"""
    num of reads: {predicitons_log["read_count"]}
    num of final multimapped reads: {predicitons_log["final_multi_gene_count"]}
    num of rescued multimapped reads: {predicitons_log["rescued_multigene_count"]}
    num of final ambiguous reads: {predicitons_log["final_ambiguous_count"]}
    """)
```

