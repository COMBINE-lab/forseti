#!/usr/bin/env python3
import os
from forseti.forseti_class import algn_to_tx_prob_parallel, ReadOriginPredictor

# NOTE:You may need to change the path to the forseti_sample_dir
forseti_sample_dir = os.getcwd()
idx_dir = os.path.join(forseti_sample_dir, 'data', 'spliceu_refs')

# load the files
reads_txome_bam_file = os.path.join(
    forseti_sample_dir, 'STAR_out', 'exonic_Aligned.toTranscriptome.out.bam')
spliceu_fasta_file = os.path.join(idx_dir, "spliceu.fa")
t2g_file = os.path.join(idx_dir, "t2g_3col.tsv")

# load the pre-built models
spline_model_file = os.path.join(
    forseti_sample_dir, "forseti", "pre_built_models", "spline_model.pkl")
mlp_model_file = os.path.join(
    forseti_sample_dir, "forseti", "pre_built_models", "mlp.pkl")

# set up the parameters
snr_min_size = 6
discount_perc = 1
polya_tail_len = 200
max_frag_len = 1000
# num_threads = 12 # default is 12
num_threads = 32

# Run the prediction pipeline
# this function will do the parallel processing of the reads
algn_to_tx_prob_results = algn_to_tx_prob_parallel(
    reads_txome_bam_file, forseti_sample_dir, spline_model_file, mlp_model_file, spliceu_fasta_file, num_threads)

# initialize the ReadOriginPredictor
rop = ReadOriginPredictor(algn_to_tx_prob_results,
                          reads_txome_bam_file, t2g_file)
prediction_dict, log = rop.run_model_pipeline()

# One example of using
for rname, forseti_read in prediction_dict.items():
    if not forseti_read.is_multi_mapped:
        prediction = forseti_read.get_predictions()
        if prediction.splicing_status == "U":
            print(forseti_read.read_name)
            print(forseti_read.is_multi_mapped)
            print(prediction.gene_id)
            print(prediction.orientation)
            print(prediction.splicing_status)
            print(prediction.max_prob)
            print(prediction.unsplice_prob)
            print(prediction.splice_prob)
            break
