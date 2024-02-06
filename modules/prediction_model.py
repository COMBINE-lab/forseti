#!/usr/bin/env python3
import os
import time
import math
import pickle
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from scipy import interpolate
from sklearn.preprocessing import OneHotEncoder
from joblib import Parallel, delayed
from forseti.forseti_class import ForsetiReads, ForsetiPredictions


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:(i + ksize)]
        kmers.append(kmer)

    return kmers


def get_algn_to_tx_prob(algn, spline, mlp, ref_seq, encoder, max_frag_len=1000,
                        polya_tail_len=200, discount_perc=1, snr_min_size=6,
                        binding_affinity_threshold=0, is_single_exon_tx=False):
    """
    This step takes read's alignment to reference transcript, spline model, mlp model and other settings, \
    and calculates the maximum of joint probabilities of arising from this reference transcript. \
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
    # frag_len_weight = 0.5
    # discount_perc = 1 if algn["reference_name"].endswith("-U") else discount_perc

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

                # we compute the binding affinity of the 30-mers that contain AAA
                # TODO: the triple A here might be confusing for users
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
    return joint_prob, best_pos


def algn_to_tx_prob_parallel(txome_bam, spline, mlp, spliceu_txome, max_frag_len, num_threads,
                             polya_tail_len=200, discount_perc=1, snr_min_size=6, binding_affinity_threshold=0):
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

    def algn_to_tx_prob(algn, spline, mlp, ref_seq, encoder, max_frag_len, polya_tail_len,
                        discount_perc, snr_min_size, binding_affinity_threshold, is_single_exon_tx):
        prob, best_pos = get_algn_to_tx_prob(algn, spline, mlp, ref_seq, encoder, max_frag_len,
                                             polya_tail_len, discount_perc, snr_min_size, binding_affinity_threshold, is_single_exon_tx)
        return (algn['query_name'], algn['reference_id'], not algn['is_reverse'], prob)

    encoder = OneHotEncoder(
        categories=[['A', 'C', 'G', 'T', 'N']] * 30, handle_unknown='ignore')
    start_time = time.time()

    algn_to_tx_prob = Parallel(n_jobs=num_threads, prefer=None)(
        delayed(algn_to_tx_prob)(
            {"is_reverse": algn.is_reverse,
             "query_length": algn.query_length,
             "cigartuples": algn.cigartuples,
             "reference_end": algn.reference_end,
             "reference_start": algn.reference_start,
             "reference_name": algn.reference_name,
             "query_name": algn.query_name,
             "reference_id": algn.reference_id
             },
            spline,
            mlp,
            spliceu_txome.get(algn.reference_name),
            encoder,
            max_frag_len,
            polya_tail_len,
            discount_perc,
            snr_min_size,
            binding_affinity_threshold,
            False
        ) for (rid, algn) in enumerate(txome_bam.fetch(until_eof=True))
    )

    elapsed_time = time.time() - start_time
    print(f"Processed the BAM file in {elapsed_time} seconds")

    return (algn_to_tx_prob)


def read_to_gene_orienataion_prob(algn_to_tx_prob, tid2g, tid2tname):
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

    return reads_predictions


def read_to_gene_prob(reads_predictions, splicing_status_min_ratio=1, antisense_min_ratio=1):
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

            is_sense = int(False)  # we might not have the final orientation
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

    return reads_gene_probs


def read_gene_prediction(reads_gene_probs, r2g=None):
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

    forseti_dict = {}
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
            forseti_dict[rname] = read_record

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
            max_index = [i for i, x in enumerate(gene_probs) if x == max_prob]
            # if we have a single best gene, then we return it
            if len(max_index) == 1:
                max_index = max_index[0]
                # save the max prediction in the ForsetiReads object
                fr_predictions = ForsetiPredictions(
                    *predictions[max_index][:4], *predictions[max_index][4])
                read_record = ForsetiReads(rname, False)
                read_record.add_prediction(fr_predictions)
                forseti_dict[rname] = read_record

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

                forseti_dict[rname] = read_record

    log = {
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

    # return read_gene_predictions, log
    return forseti_dict, log


def forseti_predictor(reads_txome_bam_file, spliceu_fasta_file, t2g_file, spline_model_file, mlp_model_file,
                      polya_tail_len, discount_perc=1, snr_min_size=6, max_frag_len=1000, num_threads=12):
    # Read in the txome bam file
    save = pysam.set_verbosity(0)
    reads_txome_bam = pysam.AlignmentFile(reads_txome_bam_file, "rb")
    pysam.set_verbosity(save)
    # Read in the t2g file

    t2g_df = pd.read_csv(t2g_file, sep="\t", header=None,
                         names=["tid", "gid", "splicing_status"])

    t2g = t2g_df.set_index("tid").to_dict()["gid"]
    tid2tname = dict(zip(range(len(reads_txome_bam.references)),
                     reads_txome_bam.references))
    tid2g = {tid: t2g[tname] for tid, tname in tid2tname.items()}

    # Read in the spline& mlp model
    with open(spline_model_file, 'rb') as file_model:
        spline = pickle.load(file_model)
    with open(mlp_model_file, 'rb') as file_model:
        mlp = pickle.load(file_model)

    # Read in the spliceu fasta file
    spliceu_txome = {record.id: record.seq for record in SeqIO.parse(
        spliceu_fasta_file, "fasta")}

    # Call the functions
    algn_to_tx_prob = algn_to_tx_prob_parallel(
        reads_txome_bam, spline, mlp, spliceu_txome, max_frag_len, num_threads,
        polya_tail_len=polya_tail_len, discount_perc=discount_perc,
        snr_min_size=snr_min_size, binding_affinity_threshold=0)

    forseti_dict, predicitons_log = read_gene_prediction(
        read_to_gene_prob(read_to_gene_orienataion_prob(algn_to_tx_prob, tid2g, tid2tname)))
    reads_txome_bam.close()

    return forseti_dict, predicitons_log
