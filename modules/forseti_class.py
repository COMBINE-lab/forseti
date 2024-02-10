#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import os
import math
import pickle
from scipy import interpolate
from Bio import SeqIO
from sklearn.preprocessing import OneHotEncoder
import time
from joblib import Parallel, delayed
import random
from pathlib import Path


def split_large_bam(large_bam_file_path, forseti_sample_dir, num_threads=12):
    """
    Split a large BAM file into smaller subfiles with a specified number of records in each.

    :param large_bam_file_path: The path to the large input BAM file to be split.
    :param forseti_sample_dir: Directory path where the split BAM files will be saved.
    :param num_threads: Number of  subfile (default is 12).
    """
    start_time = time.time()
    # Initialize counters
    batch_count = 0
    record_count = 0

    # Ensure the output directory exists
    output_dir = os.path.join(forseti_sample_dir, 'split_bam')
    os.makedirs(output_dir, exist_ok=True)

    # Open the large BAM file for reading
    large_bam_file = pysam.AlignmentFile(large_bam_file_path, "rb")

    total_algn_count = large_bam_file.count(until_eof=True)
    print(f'total_algn_count: {total_algn_count}')
    # calculate the number of records in each subfile
    records_per_file = math.ceil(
        total_algn_count/num_threads)
    print(f'records_per_file: {records_per_file}')
    # Prepare the first output file
    sub_file_path = os.path.join(output_dir, f"batch_{batch_count}.bam")
    sub_file = pysam.AlignmentFile(
        sub_file_path, "wb", template=large_bam_file, threads=num_threads)

    # Open the large BAM file for reading
    large_bam_file = pysam.AlignmentFile(large_bam_file_path, "rb")
    # Read through the input file line by line
    for read in large_bam_file.fetch(until_eof=True):
        # Write the current read to the current output file
        sub_file.write(read)
        record_count += 1

        # Check if the current output file has reached the desired size
        if record_count >= records_per_file:
            # Close the current output file
            sub_file.close()
            batch_count += 1
            record_count = 0

            # Prepare the next output file
            sub_file_path = os.path.join(
                output_dir, f"batch_{batch_count}.bam")
            sub_file = pysam.AlignmentFile(
                sub_file_path, "wb", template=large_bam_file, threads=num_threads)
    # Clean up: Close the last output file and the input file
    sub_file.close()
    large_bam_file.close()
    end_time = time.time()

    print(f"Completed. Split the BAM file into {batch_count+1} subfiles.")
    print(f"Time taken to split the bam file: {end_time - start_time} seconds")

    # top200cells: Time taken to split the bam file: 79.3 seconds


def get_algn_to_tx_prob(algn_tuple, spline, mlp, ref_seq, max_frag_len=1000,
                        polya_tail_len=200, discount_perc=1, snr_min_size=6,
                        binding_affinity_threshold=0, is_single_exon_tx=False):
    """
    This step takes read's alignment to reference transcript, spline model, mlp model and other settings, and calculates the maximum of joint probabilities of arising from this reference transcript.
    It returns the maximum joint probability and the best position of corresponding priming site.
    """

    # if this is a sense alignment w.r.t. the reference transcript
    # we first calculate the distance of the reference end site to the downstream polyA sites that are within 1000 bases of the read alignment
    """
    # x = (start, end) - polyA
    # sense reads
    #     --read-->
    # tx ---------------------------[polyA]-------------------[polyA]--> 3'
    #                    --sliding window--

    #                               --sliding window--
    #                    |          | <- final reasult window
    #
    #                --sliding window--
    #     --sliding window--
    #   [polyT]      [polyT]           <--reads--
    #                      |          | <- final reasult window
    """

    joint_prob = 0
    best_pos = 0

    algn = {"is_reverse": algn_tuple[0],
            "query_length": algn_tuple[1],
            "cigartuples": algn_tuple[2],
            "reference_end": algn_tuple[3],
            "reference_start": algn_tuple[4],
            "reference_name": algn_tuple[5],
            "query_name": algn_tuple[6],
            "reference_id": algn_tuple[7]
            }
    encoder = OneHotEncoder(
        categories=[['A', 'C', 'G', 'T', 'N']] * 30, handle_unknown='ignore')

    # test if the spline is valid
    try:
        interpolate.splev(1, spline)
    except:
        raise ValueError("The spline model is not valid.")
        spline = spline

    # test if the mlp is valid
    try:
        mlp.predict_proba(encoder.fit_transform([list("A" * 30)]))
    except:
        raise ValueError("The MLP model is not valid.")
        mlp = mlp

    def build_kmers(sequence, ksize):
        # this function takes a biopython.SeqIO.Seq and a kmer size and return a list of kmers
        kmers = []
        n_kmers = len(sequence) - ksize + 1

        for i in range(n_kmers):
            kmer = sequence[i:(i + ksize)]
            kmers.append(kmer)

        return kmers

    if not algn["is_reverse"]:
        ref_start = algn['reference_end'] - algn["query_length"]

        if algn["cigartuples"][-1][0] == 4:
            ref_start += algn["cigartuples"][-1][1]

        # if we have 30 bases downstream, we want to consider internal polyA sites
        # we use > 30 because the polyA should start one base after the read alignment
        if len(ref_seq) - ref_start > 30:
            # we want to take the 1000 30 k-mers, so it's 1030 bases
            # polyA site has to be at least 1 base away from ref start so that the ref length is positive
            downstream_30mers = build_kmers(
                ref_seq[(ref_start+1):min(ref_start + max_frag_len + 30, len(ref_seq))], 30)
            # we only compute the binding affinity of the 30-mers that contain AAAAAA
            has_enough_a = ["A" * snr_min_size in x for x in downstream_30mers]
            if has_enough_a and sum(has_enough_a) > 0:

                # first, we want to compute the fragment length probability of each 30-mer
                downstream_frag_len_prob = interpolate.splev(
                    range(
                        1,
                        len(downstream_30mers) + 1
                    ),
                    spline
                )

                # we compute the binding affinity of the 30-mers that contain AAA..
                nonzero_binding_affinity = mlp.predict_proba(encoder.fit_transform(
                    [list(x) for i, x in enumerate(downstream_30mers) if has_enough_a[i]]))[:, 1] * discount_perc
                nonzero_binding_affinity = np.array(
                    [x if x > binding_affinity_threshold else 0 for x in nonzero_binding_affinity])
                # we initialize the binding affinity as a small number
                downstream_binding_affinity = np.zeros(
                    len(downstream_30mers))  # + min(nonzero_binding_affinity)/2
                downstream_binding_affinity[np.array(
                    has_enough_a)] = nonzero_binding_affinity
                # joint_probs = downstream_binding_affinity*(1-frag_len_weight) + downstream_frag_len_prob * frag_len_weight
                joint_probs = downstream_binding_affinity * downstream_frag_len_prob

                # we want to discount the priming window probagblity to distinguish internal polyA and termial polyA
                joint_prob = joint_probs.max()
                best_pos = joint_probs.argmax() + 1

        # if we reach the end of the reference,
        # we want to consider the polyA tail

        # ----> read
        # tx ------------------------> 3'
        #           |<--------->| <- start of the last 30-mer
        #                           start pos of the last 30-mer - ref_start
        # The distance from reference start to the last kmer start
        dis_to_tx_end_30mer = len(ref_seq) - 30 + 1 - ref_start

        # if the distance to the tail is less than the max window size,
        # then we want to consider polyA tail
        # farg length range is 0-999, but the value is 1000
        if dis_to_tx_end_30mer < max_frag_len and not is_single_exon_tx:
            # we first compute the fragment length
            tail_joint_prob = interpolate.splev(
                range(
                    # one base after the start pos of the last 30-mer in the ref
                    dis_to_tx_end_30mer + 1,
                    min(max_frag_len, dis_to_tx_end_30mer + 1 + polya_tail_len) + 1
                ),
                spline
            )  # * frag_len_weight
            # tail_binding_affinity = np.ones(min(max_frag_len - 30 + 1 - dist_to_tail, 200))
            tail_30mers = build_kmers(
                ref_seq[-(30 - 1):] + 'A' * min(max_frag_len - dis_to_tx_end_30mer, 29), 30)
            # for others, we use the binding affinity of all A
            # all_a_prob = mlp.predict_proba(encoder.fit_transform([list('A'*30)]))[:,1]
            all_a_prob = 1
            tail_binding_affinity = mlp.predict_proba(
                encoder.fit_transform([list(x) for x in tail_30mers]))[:, 1]
            tail_binding_affinity = np.array(
                [x if x > binding_affinity_threshold else 0 for x in tail_binding_affinity])

            # tail_binding_affinity = np.array([1 if x > 0.7 else 0 for x in tail_binding_affinity])

            tail_joint_prob[:len(tail_30mers)] = tail_joint_prob[:len(
                tail_30mers)] * tail_binding_affinity  # + tail_binding_affinity * (1-frag_len_weight)
            tail_joint_prob[len(tail_30mers):] = tail_joint_prob[len(
                tail_30mers):] * all_a_prob  # + all_a_prob #* (1-frag_len_weight)

            max_tail_joint_prob = tail_joint_prob.max()
            # if we get the maximum probability from the tail, we set the best position to be infinity
            if max_tail_joint_prob > joint_prob:
                joint_prob = max_tail_joint_prob
                best_pos = float('inf')

    else:
        # now we have a antisense alignment
        """
        antisense reads
            read                                                                       <----------
              tx   ========================================================================================================> 3'
                                                      [polyT)
                                       ------------------------------------
        """

        ref_end = algn["reference_start"] + algn['query_length']

        if algn["cigartuples"][0][0] == 4:
            ref_end += algn['cigartuples'][0][1]

        # if we have 30 bases upstream, we want to consider internal polyT sites
        if ref_end > 30:
            # we take the upstream 1000 30 k-mers
            upstream_30mers = build_kmers(
                ref_seq[max(ref_end - max_frag_len - 30, 0):ref_end - 1].reverse_complement(), 30)

            has_enough_a = ["A" * snr_min_size in x for x in upstream_30mers]
            if has_enough_a and sum(has_enough_a) > 0:
                # first, we want to compute the fragment length probability of each 30-mer
                upstream_frag_len_prob = interpolate.splev(
                    range(
                        1,
                        1 + len(upstream_30mers)
                    ),
                    spline
                )

                # we initialize the binding affinity as zero
                # we compute the binding affinity of the 30-mers that contain AAA
                nonzero_binding_affinity = mlp.predict_proba(encoder.fit_transform(
                    [list(x) for i, x in enumerate(upstream_30mers) if has_enough_a[i]]))[:, 1] * discount_perc
                nonzero_binding_affinity = np.array(
                    [x if x > binding_affinity_threshold else 0 for x in nonzero_binding_affinity])

                # + min(nonzero_binding_affinity)/2
                upstream_binding_affinity = np.zeros(len(upstream_30mers))
                upstream_binding_affinity[np.array(
                    has_enough_a)] = nonzero_binding_affinity

                # joint_probs = upstream_binding_affinity * (1-frag_len_weight) + upstream_frag_len_prob * frag_len_weight
                joint_probs = upstream_binding_affinity * upstream_frag_len_prob

                # we want to discount the priming window probagblity to distinguish internal polyA and termial polyA
                joint_prob = joint_probs.max()
                best_pos = joint_probs.argmax() + 1

    # Then we return the maximum probability
    return algn['query_name'], algn['reference_id'], algn['is_reverse'], joint_prob, best_pos


def algn_to_tx_prob_parallel(reads_txome_bam_file, forseti_sample_dir,  spline_model_file, mlp_model_file, spliceu_fasta_file, num_threads=12, max_frag_len=1000, polya_tail_len=200, discount_perc=1, snr_min_size=6, binding_affinity_threshold=0):
    """
    This function is used to parallelize the `algn_to_tx_prob` function
    It takes
    1. a pysam.AlignmentFile object from the transcriptome BAM file
    2. the polyA and polyT sites on each transcript
    3. the spline model
    4. the MLP model
    5. the spliceu reference sequence
    6. the polyA minimum length
    7. the maximum fragment length
    8. the number of threads to use

    It returns a list of tuples, where each tuple is \\
    (read name, reference id, is sense, maximum joint probability of arising from the reference).\\
    Each alignment in the BAM file will be converted to a tuple after computing the probability of arising from the reference.
    """
    # Use joblib's Parallel to parallelize predictions
    # cannot use prefer="threads" as it will move the variable into the thread

    def algn_to_tx_prob(txome_bam_path, spline, mlp, spliceu_txome, max_frag_len, polya_tail_len, discount_perc, snr_min_size, binding_affinity_threshold, is_single_exon_tx):
        # Read in the txome bam file
        save = pysam.set_verbosity(0)
        txome_bam = pysam.AlignmentFile(txome_bam_path, "rb")
        pysam.set_verbosity(save)

        algn_list = [(algn.is_reverse,
                      algn.query_length,
                      algn.cigartuples,
                      algn.reference_end,
                      algn.reference_start,
                      algn.reference_name,
                      algn.query_name,
                      algn.reference_id)
                     for algn in txome_bam.fetch(until_eof=True)]
        # some reads may take very long time to process, shuffle the list to avoid the all alignments to those reads are grouped together.
        random.seed(42)
        random.shuffle(algn_list)
        algn_to_tx_prob_results = []
        for algn_tuple in algn_list:
            ref_seq = spliceu_txome.get(algn_tuple[5])
            algn_query_name, algn_reference_id, algn_is_reverse, prob, best_pos = get_algn_to_tx_prob(
                algn_tuple, spline, mlp, ref_seq, max_frag_len, polya_tail_len, discount_perc, snr_min_size, binding_affinity_threshold, is_single_exon_tx)
            algn_to_tx_prob_results.append(
                (algn_query_name, algn_reference_id, not algn_is_reverse, prob))
        return algn_to_tx_prob_results

    # split large bam file into smaller subfiles
    # new files will be saved in the 'split_bam' folder under the 'forseti_sample_dir'
    split_large_bam(reads_txome_bam_file, forseti_sample_dir, num_threads)
    split_bam_folder_path = os.path.join(forseti_sample_dir, 'split_bam')
    sub_bam_list = [file for file in Path(split_bam_folder_path).iterdir()
                    if file.is_file()]

    # we first load the two built-in models
    try:
        with open(spline_model_file, 'rb') as file_model:
            spline = pickle.load(file_model)
    except:
        raise ValueError("The spline model is not valid.")
    try:
        with open(mlp_model_file, 'rb') as file_model:
            mlp = pickle.load(file_model)
    except:
        raise ValueError("The MLP model is not valid.")
    # prepare the spliceu txome map
    spliceu_txome = {record.id: record.seq for record in SeqIO.parse(
        spliceu_fasta_file, "fasta")}

    start_time = time.time()
    all_bam_result = Parallel(n_jobs=num_threads, prefer=None)(delayed(algn_to_tx_prob)(txome_bam_path, spline, mlp, spliceu_txome,
                                                                                        max_frag_len, polya_tail_len, discount_perc, snr_min_size, binding_affinity_threshold, False) for txome_bam_path in sub_bam_list)

    elapsed_time = time.time() - start_time
    print(f"Processed the BAM file(s) in {elapsed_time} seconds")
    algn_to_tx_prob = [
        algn_result for sub_bam_result in all_bam_result for algn_result in sub_bam_result]
    return algn_to_tx_prob


class ReadOriginPredictor:
    # TODO: modify the description

    """
    define a function that takes a read alignment, and the polyA intervals of the reference, 
    and return the maximum probability of the read arising from the reference
    - Input: 
        1. a pysam.AlignedSegment
        2. a list of tuples, where each tuple is a polyA interval on the reference
    - Output:
        1. the maximum probability of the read arising from the reference

    ## details
    In this model, we process a transcriptomic read alignment (read <-> tx, plus some affiliate information) 
    to find the probability that the read arises from the reference transcript.
    The idea is, 
    1. we have a cDNA fragment length distribution.
    2. we know where are the potential priming sites, the polyA sites on the reference transcript. 
    3. We know that the interval spanning the priming site and the read mapping site represene the cDNA fragment length.

    Therefore, giving a read alignment, we can calculate the size of the interval between the read alignment and all possible polyA sites. 
    Then, we apply the cDNA fragment length distribution to 
    convert the interval size to the probability of the read arising from the reference transcript.

    The difficulty is that although we know there should be a polyA in the priming window, we do not know the exact location of the polyA.
     Therefore, we need to consider all possible priming window containing the polyA site and select the one with the maximum probability. 

    """

    def __init__(self, algn_to_tx_prob_results, reads_txome_bam_file, t2g_path):
        """
        Initialize the ReadOriginPredictor class

        Args:
            transcriptome_bam_path (str): A str to the transcriptome BAM file
            spliceu_path (str): path to the spliceu reference sequence 
            spline (scipy.interpolate): the spline model for the fragment length distribution
            mlp (sklearn.neural_network.MLPClassifier): the MLP model for the priming site prediction
            t2g (dict): the transcript to gene mapping
            max_frag_len (int): the maximum fragment length
            discount_perc (float): the discount percentage for internal polyA sites
            polya_tail_len (int): the length of the polyA tail
        """
        self.algn_to_tx_prob_results = algn_to_tx_prob_results

        save = pysam.set_verbosity(0)
        reads_txome_bam = pysam.AlignmentFile(reads_txome_bam_file, "rb")
        pysam.set_verbosity(save)
        # we get tid2tname
        self.tid2tname = dict(
            zip(range(len(reads_txome_bam.references)), reads_txome_bam.references))

        # we read in the t2g file
        t2g_df = pd.read_csv(t2g_path, sep="\t", header=None,
                             names=["tid", "gid", "splicing_status"])
        self.t2g = t2g_df.set_index("tid").to_dict()["gid"]
        self.tid2g = {tid: self.t2g[tname]
                      for tid, tname in self.tid2tname.items()}
        # close the BAM file
        reads_txome_bam.close()

        self.read_to_gene_orientation_results = None
        self.read_to_gene_results = None
        self.forseti_prediction_dict = None
        self.forseti_prediction_log = None

    def read_to_gene_orienataion_prob(self):
        """
        This function takes the alignment to transctipt probability returned from `algn_to_tx_prob_parallel` and convert it to read to gene probability.
        The input object is a list of tuples, where each tuple is \\
        (read name, reference id, is sense, probability of arising from the reference). \\
        Each read in the BAM file will be converted to a tuple after summarizing all its alignments.
        This function outputs a dictionary. The dictionary contains the read to gene probability for each reference gene.
        Specifically, for each read, to each reference gene, for both sense and antisense alignments, \
            and for both splice and unsplice, we keep the maximum probability, the sum probability, and the number of appearances.

        Input:  1. a list of tuples:
                    (read name, reference id, is sense, maximum joint probability of arising from the reference).
                2. tid2g: tanscript id to gene mapping
                3. tid2tname : transcript id to transcript name mapping

        Output: a dictionary of
        {
            read_name: {
                gid: [
                    [antisense_spliced_prob_sum, antisense_unspliced_prob_sum, antisense_max_spliced_prob, antisense_max_unspliced_prob, num_unspliced, num_spliced],
                    [sense_spliced_prob_sum, sense_unspliced_prob_sum, sense_max_spliced_prob, sense_max_unspliced_prob, num_unspliced, num_spliced]
                ]
            }
        } 
        """
        # keys: read name
        # values: {gid: [2*[unspliced_prob_sum, spliced_prob_sum, max_unspliced_prob, max_spliced_prob, num_unspliced, num_spliced]]}
        # where the first list is for antisense alignments and the second list is for sense alignments

        algn_to_tx_prob = self.algn_to_tx_prob_results
        tid2g = self.tid2g
        tid2tname = self.tid2tname

        # we check if we have a valid algn_to_tx_prob
        if algn_to_tx_prob is None:
            raise ValueError("Please run algn_to_tx_prob_parallel first.")
        reads_predictions = {it[0]: {} for it in algn_to_tx_prob}
        for rname, tid, is_sense, prob in algn_to_tx_prob:
            gid = tid2g[tid]

            # if this is a new reference gene
            if gid not in reads_predictions[rname]:
                # if we see the same best probability, we add it in
                reads_predictions[rname][gid] = [
                    [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

            # if the probability is zero, we skip it
            if prob == 0:
                continue

            # we define the index of the list we want to modify
            is_sense = int(is_sense)
            is_spliced = int(not tid2tname[tid].endswith("-U"))
            # we update the sum probability
            reads_predictions[rname][gid][is_sense][is_spliced] += prob
            # we update the number of appearances
            reads_predictions[rname][gid][is_sense][is_spliced+4] += 1
            # we update the max probability
            if prob > reads_predictions[rname][gid][is_sense][is_spliced+2]:
                reads_predictions[rname][gid][is_sense][is_spliced+2] = prob

        self.read_to_gene_orientation_results = reads_predictions
        return reads_predictions

    def read_to_gene_prob(self, splicing_status_min_ratio=1, antisense_min_ratio=1):
        """
        This step takes the read to gene probability returned from `read_to_gene_orienataion_prob` and convert it to `read to gene probability`.
        Input: a dictionary of
        {
            read_name: {
                gid: [
                    [antisense_spliced_prob_sum, antisense_unspliced_prob_sum, antisense_max_spliced_prob, antisense_max_unspliced_prob, num_unspliced, num_spliced],
                    [sense_spliced_prob_sum, sense_unspliced_prob_sum, sense_max_spliced_prob, sense_max_unspliced_prob, num_unspliced, num_spliced]
                ]
            }
        } 
        Output: a dictionary of
        {
            read_name: [[
                orientation,
                splicing_status,
            ]]
        }

        """
        def get_splicing_status(prob_list, min_ratio):
            """
            In this function, we assign splicing status based on the max s/u probability and require high max s/u prob ratio.
            """

            # [sum, max, num]
            is_spliced = 0
            splicing_status = "A"
            s_prob = 0.5
            u_prob = 0.5

            if sum(prob_list) == 0:
                return (splicing_status, 0, 0, (0.5, 0.5))

            # get mean from sum (the first twos), and count (the third twos)
            mean_u_prob = prob_list[0]/prob_list[4] if prob_list[4] != 0 else 0
            mean_s_prob = prob_list[1]/prob_list[5] if prob_list[5] != 0 else 0
            max_u_prob = prob_list[2]
            max_s_prob = prob_list[3]

            # if sum probs are the same, we use the winner of the max probabilities
            if max_u_prob != max_s_prob:
                is_spliced = int(max_s_prob > max_u_prob)
                splicing_status = "S" if is_spliced == 1 else "U"
                denominator = max_s_prob + max_u_prob
                s_prob = max_s_prob/denominator
                u_prob = max_u_prob/denominator

            # # if we have a winner, we use it
            # elif mean_u_prob != mean_s_prob:
            #     is_spliced = int(mean_s_prob > mean_u_prob)
            #     splicing_status = "S" if is_spliced == 1 else "U"
            #     denominator = mean_s_prob + mean_u_prob
            #     s_prob = mean_s_prob/denominator
            #     u_prob = mean_u_prob/denominator

            mean_prob = mean_s_prob if is_spliced == 1 else mean_u_prob
            max_prob = prob_list[is_spliced+2]

            # # finally, we want to change the splicing status back to ambiguous if the ratio is too small
            # if mean_s_prob == 0 or mean_u_prob == 0:
            #     mean_prob_ratio = min_ratio
            # else:
            #     mean_prob_ratio = mean_s_prob/mean_u_prob if is_spliced == 1 else mean_u_prob/mean_s_prob

            if max_s_prob == 0 or max_u_prob == 0:
                max_prob_ratio = min_ratio
            else:
                max_prob_ratio = max_s_prob/max_u_prob if is_spliced == 1 else max_u_prob/max_s_prob

            # if mean_prob_ratio < min_ratio and max_prob_ratio < min_ratio:
            if max_prob_ratio < min_ratio:
                splicing_status = "A"

            return (splicing_status, mean_prob, max_prob, (u_prob, s_prob))

        reads_predictions = self.read_to_gene_orientation_results
        if reads_predictions is None:
            raise ValueError("Please run read_to_gene_orienataion_prob first.")

        reads_gene_probs = {rname: [] for rname in reads_predictions}
        for rname, predictions in reads_predictions.items():
            # predictions = {gid: [2*[unspliced_prob_sum, spliced_prob_sum, max_unspliced_prob, max_spliced_prob, num_unspliced, num_spliced]]}
            for gid, probs in predictions.items():
                # we first process antisense alignments
                is_sense = int(False)
                antisense_splicing_status, antisense_mean_prob, antisense_max_prob, antisense_prediction_prob = get_splicing_status(
                    probs[is_sense], min_ratio=splicing_status_min_ratio)

                # we then process sense alignments
                is_sense = int(True)
                sense_splicing_status, sense_mean_prob, sense_max_prob, sense_prediction_prob = get_splicing_status(
                    probs[is_sense], min_ratio=splicing_status_min_ratio)

                # we might not have the final orientation
                is_sense = int(False)
                splicing_status = "A"  # initialize as ambiguous
                orientation = "*"
                final_prediction_prob = (0.5, 0.5)

                # if we have a winner, we use it
                # If we are totally ambiguous
                if antisense_mean_prob == sense_mean_prob and antisense_max_prob == sense_max_prob:
                    reads_gene_probs[rname].append(
                        (gid, orientation, splicing_status, antisense_max_prob, final_prediction_prob))
                    continue
                # if not, we use the winner of the max probabilities
                if antisense_max_prob != sense_max_prob:
                    is_sense = int(sense_max_prob >
                                   antisense_max_prob/antisense_min_ratio)
                    splicing_status = sense_splicing_status if is_sense == 1 else antisense_splicing_status
                    final_prediction_prob = sense_prediction_prob if is_sense == 1 else antisense_prediction_prob
                    orientation = "+" if is_sense == 1 else "-"
                elif antisense_mean_prob != sense_mean_prob:
                    # we consider it as antisense only if its antisense prob is antisense_min_ratio(i.e. 1.2) times larger than sense
                    is_sense = int(sense_mean_prob >
                                   antisense_mean_prob/antisense_min_ratio)
                    splicing_status = sense_splicing_status if is_sense == 1 else antisense_splicing_status
                    final_prediction_prob = sense_prediction_prob if is_sense == 1 else antisense_prediction_prob
                    orientation = "+" if is_sense == 1 else "-"

                max_prob = sense_max_prob if is_sense == 1 else antisense_max_prob

                # Now we push the result to the dictionary
                reads_gene_probs[rname].append(
                    (gid, orientation, splicing_status, max_prob, final_prediction_prob))

        self.read_to_gene_results = reads_gene_probs
        return reads_gene_probs

    def read_gene_prediction(self, r2g=None):
        """
        This function takes the read to its reference genes' probability returned from `read_to_gene_prob` and convert it to read to gene prediction. \
            We will return a list of references that the read is most likely to arise from. If there is a tie, we will return all the references.
        Input: a dictionary of
        {
            read_name: [
                (gid, orientation, splicing_status, max_prob, (u_prob, s_prob))
            ]
        }

        Output: a dictionary of
        {
            read_name: [
                (gid, orientation, splicing_status, max_prob, (u_prob, s_prob))
            ]
        }
        """
        read_count = 0
        unigene_count = 0
        unigene_spliced_reads_count = 0
        unigene_unspliced_reads_count = 0
        unigene_ambiguous_reads_count = 0
        unigene_wronggene_count = 0
        multigene_count = 0
        rescued_multigene_count = 0
        rescued_wrong_gene_count = 0
        rescued_spliced_reads_count = 0
        rescued_unspliced_reads_count = 0
        rescued_ambiguous_reads_count = 0
        final_multi_gene_count = 0

        reads_gene_probs = self.read_to_gene_results
        forseti_prediction_dict = {}
        for rname, predictions in reads_gene_probs.items():
            # predictions : [(gid, orientation, splicing_status, max_prob, (unsplice_prob, splice_prob))]
            read_count += 1
            # if we have a single gene, we just return it
            if len(predictions) == 1:

                # save the prediction in the ForsetiReads object
                fr_predictions = ForsetiPredictions(
                    *predictions[0][:4], *predictions[0][4])
                read_record = ForsetiReads(rname, False)
                read_record.add_prediction(fr_predictions)
                forseti_prediction_dict[rname] = read_record

                unigene_count += 1
                # if the best gene is not the true gene(applicable with simulation data/paired-end data\
                # with known origin/inferred origin), we count it as wrong
                if r2g is not None and predictions[0][0] != r2g[rname]:
                    unigene_wronggene_count += 1
                elif predictions[0][2] == "S":
                    unigene_spliced_reads_count += 1
                elif predictions[0][2] == "U":
                    unigene_unspliced_reads_count += 1
                else:
                    unigene_ambiguous_reads_count += 1
            # now we need to do a prediction
            else:
                multigene_count += 1
                # get genes' probability
                gene_probs = [x[3] for x in predictions]
                gene_oris = [x[1] for x in predictions]
                max_prob = max(gene_probs)
                # get the best genes
                max_index = [i for i, x in enumerate(
                    gene_probs) if x == max_prob]
                # if we have a single best gene, then we return it
                if len(max_index) == 1:
                    max_index = max_index[0]
                    # save the max prediction in the ForsetiReads object
                    fr_predictions = ForsetiPredictions(
                        *predictions[max_index][:4], *predictions[max_index][4])
                    read_record = ForsetiReads(rname, False)
                    read_record.add_prediction(fr_predictions)
                    forseti_prediction_dict[rname] = read_record

                    rescued_multigene_count += 1
                    # if the best gene is not the true gene, we count it as wrong
                    if r2g is not None and predictions[max_index][0] != r2g[rname]:
                        rescued_wrong_gene_count += 1
                    elif predictions[max_index][2] == "S":
                        rescued_spliced_reads_count += 1
                    elif predictions[max_index][2] == "U":
                        rescued_unspliced_reads_count += 1
                    else:
                        rescued_ambiguous_reads_count += 1
                # if we have multiple best genes, we return all of them
                else:
                    final_multi_gene_count += 1
                    read_record = ForsetiReads(rname, True)
                    # save the all predictions into ForsetiReads object
                    for i in max_index:
                        read_record.add_prediction(
                            ForsetiPredictions(*predictions[i][:4], *predictions[i][4]))

                    forseti_prediction_dict[rname] = read_record

        forseti_prediction_log = {
            "read_count": read_count,
            "unigene_count": unigene_count,
            "unigene_wronggene_count": unigene_wronggene_count,
            "unigene_spliced_reads_count": unigene_spliced_reads_count,
            "unigene_unspliced_reads_count": unigene_unspliced_reads_count,
            "unigene_ambiguous_reads_count": unigene_ambiguous_reads_count,
            "multigene_count": multigene_count,
            "rescued_wrong_gene_count": rescued_wrong_gene_count,
            "rescued_multigene_count": rescued_multigene_count,
            "rescued_spliced_reads_count": rescued_spliced_reads_count,
            "rescued_unspliced_reads_count": rescued_unspliced_reads_count,
            "rescued_ambiguous_reads_count": rescued_ambiguous_reads_count,
            "final_multi_gene_count": final_multi_gene_count,
            "final_spliced_count": unigene_spliced_reads_count + rescued_spliced_reads_count,
            "final_unspliced_count": unigene_unspliced_reads_count + rescued_unspliced_reads_count,
            "final_ambiguous_count": unigene_ambiguous_reads_count + rescued_ambiguous_reads_count
        }

        self.forseti_prediction_dict = forseti_prediction_dict
        self.forseti_prediction_log = forseti_prediction_log

        return forseti_prediction_dict, forseti_prediction_log

    def run_model_pipeline(self):
        print('start: run_model_pipeline')

        self.read_to_gene_orienataion_prob()
        self.read_to_gene_prob()
        forseti_prediction_dict, forseti_prediction_log = self.read_gene_prediction()
        return forseti_prediction_dict, forseti_prediction_log


class ForsetiReads:
    def __init__(self, read_name, is_multi_mapped, forseti_predictions=None):
        self.read_name = read_name
        self.is_multi_mapped = is_multi_mapped
        # Initialize predictions as an empty list if None is provided
        self.predictions = forseti_predictions if forseti_predictions is not None else []

    def add_prediction(self, prediction):
        """Add a ForsetiPredictions object to the predictions list."""
        self.predictions.append(prediction)

    def get_predictions(self):
        """Return the list of ForsetiPredictions objects."""
        if self.is_multi_mapped:
            return self.predictions
        else:
            return self.predictions[0]

    # If needed, you can add more methods to remove or update predictions


class ForsetiPredictions:
    def __init__(self, gene_id, orientation, splicing_status, max_prob, unsplice_prob, splice_prob):
        self.gene_id = gene_id
        self.orientation = orientation
        self.splicing_status = splicing_status
        self.max_prob = max_prob
        self.unsplice_prob = unsplice_prob
        self.splice_prob = splice_prob
