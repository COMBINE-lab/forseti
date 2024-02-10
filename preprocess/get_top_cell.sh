# Please NOTE: You may need to change the path to the FST_SAMPLE_DIR!!
export FST_SAMPLE_DIR=$PWD
export FASTQ_DIR="$FST_SAMPLE_DIR/data/pbmc_1k_v3_fastqs"
export REF_DIR="$FST_SAMPLE_DIR/data/refdata-gex-GRCh38-2020-A"
export IDX_DIR="$FST_SAMPLE_DIR/data/spliceu_refs"
export STAR_IDX_DIR="$FST_SAMPLE_DIR/data/star_spliceu_refs"

export STAR_out=$FST_SAMPLE_DIR/STAR_out
mkdir -p $STAR_out

STAR --runThreadN 12 \
    --readFilesIn $FASTQ_DIR/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,$FASTQ_DIR/pbmc_1k_v3_S1_L002_R2_001.fastq.gz $FASTQ_DIR/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,$FASTQ_DIR/pbmc_1k_v3_S1_L002_R1_001.fastq.gz\
    --readFilesCommand zcat \
    --genomeDir $STAR_IDX_DIR \
    --outFileNamePrefix $FST_SAMPLE_DIR/STAR_out/ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --soloType CB_UMI_Simple \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist $FST_SAMPLE_DIR/data/whitelist/3M-february-2018.txt \
    --soloFeatures GeneFull \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --limitIObufferSize 50000000 50000000 



# Optional: For this tutorial, we extract the top cells with at least 2000 genes for this tutorial. The filtering is done based on the gene count matrix provided on 10X website. The top cells are saved in top500cells.txt.

# Download the gene count matrix folder
wget -qO- https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz | tar xzf - -C $STAR_out

# Extract the CellBarcode for top cells(have > 2000 genes)
Python extract_top_cells.py $STAR_out/filtered_feature_bc_matrix $STAR_out
# set the path to the toGenome bam file. (only the toGenome bam file has the corrected cell barcode, save in tag 'CB'; toTxome bam file has the raw cell barcode, save in tag 'CR')
BAM_FILE=$STAR_out/Aligned.sortedByCoord.out.bam
# Save the header lines
samtools view -H $BAM_FILE > $STAR_out/SAM_header
# Filter alignments using filter.txt. 
samtools view $BAM_FILE | grep -F -f $STAR_out/top200cells.txt > $STAR_out/filtered_SAM_body

# Combine header and body
cat $STAR_out/SAM_header $STAR_out/filtered_SAM_body > $STAR_out/CB_filtered.sam

# Convert filtered.sam to BAM format
samtools view -b $STAR_out/CB_filtered.sam > $STAR_out/CB_filtered.bam

# The following commands are the same as the standard verion in tutorial
samtools view -d NH:1 $STAR_out/CB_filtered.bam -o $STAR_out/filtered_Aligned.sortedByCoord.out.bam

# extract the list of unique map reads
samtools view $STAR_out/filtered_Aligned.sortedByCoord.out.bam | cut -f1 > $STAR_out/unimap_read_names.txt
# extract the toTx record based on the unique maped read list
samtools view \
    -N $STAR_out/unimap_read_names.txt \
    -o $STAR_out/filtered_Aligned.toTranscriptome.out.bam \
    $STAR_out/Aligned.toTranscriptome.out.bam

# find the reads 100% mapped in a exonic region (i.e. ambiguous reads)
bedtools intersect \
-a $STAR_out/filtered_Aligned.sortedByCoord.out.bam \
-b $IDX_DIR/exon_by_tx.bed \
-f 1.0 \
-wo \
-bed \
> $STAR_out/exonic_reads.bed

# extract the list of exonic reads
cut -f 4 $STAR_out/exonic_reads.bed > $STAR_out/exonic_read_name.txt
# filter the exonic reads in toTranscriptome.bam file
samtools view -N $STAR_out/exonic_read_name.txt -o $STAR_out/exonic_Aligned.toTranscriptome.out.bam $STAR_out/filtered_Aligned.toTranscriptome.out.bam
