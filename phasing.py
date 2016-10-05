#!/usr/bin/env python3

import scipy.sparse
import numpy as np
import scipy.stats
import sys
import logging


def matrix_to_str(np_matrix, n_rows=np.inf, n_cols=np.inf):
    '''
    Convert a subsection of a np matrix to a nice table.
    '''
    if n_cols > np_matrix.shape[1]:
        n_cols = np_matrix.shape[1]
    if n_rows > np_matrix.shape[0]:
        n_rows = np_matrix.shape[0]

    out_table = ''
    for index, v in np.ndenumerate(np_matrix[:n_rows, :n_cols].todense()):
        if index[1] == 0 and index[0] > 0:
            out_table += '\n'
        elif index[1] > 0:
            out_table += '\t'
        out_table += str(v)
    return out_table


def phase_mosaic_var(variant_id, variant_barcodes, barcodes):
    '''
    Performs phasing of the mosaic variant to nearby germline variants.
    Returns the number of barcodes used in the phasing, \
    the p-value of the variant as a mosaic and \
    the p-value of the variant being a false-positive.
    '''
    n_variants = sum([x["is_germline"] for x in variant_barcodes.values()])
    n_barcodes = len(barcodes)
    (variant_matrix, variant_index, barcode_index) = \
        construct_germline_barcode_matrix(
        variant_barcodes, n_variants, n_barcodes)
    logging.debug("Analyzed {} variants and {} barcodes".format(
            len(variant_index), len(barcode_index)))
    logging.debug("Table is:\n{}".format(
            matrix_to_str(variant_matrix, 10000, 50)))
    (haplotypes, confidence, n_seen) = get_haplotypes(variant_matrix)
    skip_barcode_indices = set()
    barcode_not_seen = set()
    barcode_discordant = set()
    for barcode_iter in range(0, n_barcodes):
        if n_seen[barcode_iter] == 0:
            barcode_not_seen.add(barcode_iter)
        elif confidence[barcode_iter] / n_seen[barcode_iter] < 0.5:
            barcode_discordant.add(barcode_iter)
    logging.debug("Removed {} barcodes due to discordance".format(
            len(barcode_discordant)))
    logging.debug("Removed {} barcodes due to not seen".format(
            len(barcode_not_seen)))
    skip_barcode_indices = barcode_not_seen | barcode_discordant
    return determine_mosaicism(haplotypes, skip_barcode_indices, barcodes,
                               barcode_index, variant_id, variant_barcodes)


def construct_germline_barcode_matrix(variant_barcodes, n_variants,
                                      n_barcodes, dtype=np.dtype("int32")):
    '''
    Construct a sparse matrix of the variants \
    supported by the barcode information.
    '''
    variant_barcode_matrix = scipy.sparse.lil_matrix(
            (n_variants, n_barcodes), dtype=np.dtype("int32"))
    variant_index = {}
    barcode_index = {}
    variant_index_counter = 0
    barcode_index_counter = 0
    for variant, properties in variant_barcodes.items():
        if not properties["is_germline"]:  # Skip non-het germline vars
            continue
        variant_index[variant] = variant_index_counter
        variant_index_counter += 1
        for allele, barcode_list in properties.items():
            if not (type(allele) is int and type(barcode_list) is list):
                continue
            for barcode in barcode_list:
                if barcode not in barcode_index:
                    barcode_index[barcode] = barcode_index_counter
                    barcode_index_counter += 1
                variant_barcode_matrix[variant_index[variant],
                                       barcode_index[barcode]] = allele + 1
    return (variant_barcode_matrix.tocsr(), variant_index, barcode_index)


def get_haplotypes(variant_matrix):
    '''
    Determine the haplotypes of the barcodes.
    Returns the haplotype (0 or 1) the number of variants \
    on the barcode and the haplotype confidence \
    (concordant variants - discordant variants).
    '''
    n_barcodes = variant_matrix.shape[1]
    haplotypes = [0] * n_barcodes
    confidence = [0] * n_barcodes  # n_concordant - n_discordant
    n_seen = [0] * n_barcodes  # Total number of times barcode was seen
    for variant_iter in range(0, variant_matrix.shape[0]):

        # Flip the variant's alleles if the data is #
        #  more consistant with the other haplotype #
        n_discordant = 0
        n_counted = 0
        is_flipped = False
        for barcode_iter in range(0, n_barcodes):
            current = variant_matrix[variant_iter, barcode_iter]
            if current:
                n_counted += 1
                if current - 1 != haplotypes[barcode_iter]:
                    n_discordant += 1
        if n_discordant > n_counted / 2:
            is_flipped = True

        # Add the variant information to the haplotypes #
        for barcode_iter in range(0, n_barcodes):
            current = variant_matrix[variant_iter, barcode_iter]
            if not current:  # barcode not represented at position
                continue
            n_seen[barcode_iter] += 1
            if ((not is_flipped and haplotypes[barcode_iter] == current - 1) or
                    (is_flipped and haplotypes[barcode_iter] != current - 1)):
                confidence[barcode_iter] += 1
            else:
                confidence[barcode_iter] -= 1
                if confidence[barcode_iter] < 0:
                    confidence[barcode_iter] = abs(confidence[barcode_iter])
                    haplotypes[barcode_iter] = current - 1
    return (haplotypes, confidence, n_seen)


def determine_mosaicism(haplotypes, skip_barcode_indices, barcodes2,
                        barcode_index, variant_id, variant_barcodes,
                        sequencing_error_rate=0.05):
    '''
    Given the barcodes that the haplotypes are a part of, \
    identify mosaic variants.
    Returns the number of barcodes used in the phasing, \
    the probability of the variant as a mosaic and the probability \
    of the variant begin a false-positive
    '''

    '''
    Phasing evidence is stored in a 2x2 table. \
    Axis 0 defines the status of the barcode, \
    axis 1 defines the status of the mosaic variant.
    '''
    phasing_evidence = np.array([[0, 0], [0, 0]])
    variant_properties = variant_barcodes[variant_id]
    noob_dp = 0
    for allele, barcodes in variant_properties.items():
        if not (type(allele) is int and type(barcodes) is list):
            continue
        for barcode in barcodes:
            # Skip barcodes not in germline variants #
            if barcode not in barcode_index:
                continue
            cur_index = barcode_index[barcode]
            if cur_index in skip_barcode_indices:
                continue
            haplotype = haplotypes[cur_index]
            phasing_evidence[haplotype, allele] += 1
            for var in barcodes2[barcode]:
                if var[:-2] in variant_barcodes:
                    if variant_barcodes[var[:-2]]["is_germline"]:
                        noob_dp += 1
                        found
    logging.debug("DP should be {}".format(noob_dp))
    logging.debug("Constructed phasing matrix {}".format(
            str(phasing_evidence)))

    depth = sum(sum(phasing_evidence))
    mosaic_haplotype = 0  # The haplotype with most mosaic variant barcodes
    if phasing_evidence[1, 1] > phasing_evidence[0, 1]:
        mosaic_haplotype = 1
    # Skip sites with no haplotype informaive reads
    #  or data on the mosaic haplotype
    prob_mosaic, prob_false_positive = 0, 0
    if (depth == 0 or (phasing_evidence[1, 1] == 0 and
                       phasing_evidence[0, 1] == 0)):
        prob_mosaic, prob_false_positive = 0, 0
    else:
        if phasing_evidence[mosaic_haplotype, 1] > 2:
            if phasing_evidence[mosaic_haplotype, 0] > 2:
                prob_mosaic = 0.98
            if phasing_evidence[abs(mosaic_haplotype - 1), 1] > 2:
                prob_false_positive = 0.98

    return (depth, prob_mosaic, prob_false_positive)
