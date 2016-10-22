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
    for index, v in np.ndenumerate(np_matrix[:n_rows, :n_cols]):
        if index[1] == 0 and index[0] > 0:
            out_table += '\n'
        elif index[1] > 0:
            out_table += '\t'
        out_table += str(v)
    return out_table


def phase_mosaic_var(variant_id, variant_barcodes):
    '''
    Performs phasing of the mosaic variant to nearby germline variants.
    Returns the number of barcodes used in the phasing, \
    the p-value of the variant as a mosaic and \
    the p-value of the variant being a false-positive.
    '''
    (variant_matrix, variant_index, barcode_index) = \
        construct_germline_barcode_matrix(variant_barcodes)

    logging.debug("Analyzed {} variants and {} barcodes".format(
            len(variant_index), len(barcode_index)))
#    logging.debug("Table is:\n{}".format(
#            matrix_to_str(variant_matrix, 10000, 50)))
    (barcode_haplotypes, barcode_concord, variant_haplotypes, variant_concord, n_seen) = get_haplotypes_DP(variant_matrix)
    #(haplotypes, confidence, n_seen) = get_haplotypes_DP(variant_matrix)
    
    skip_barcode_indices = set()
    barcode_not_seen = set()
    barcode_discordant = set()
    for barcode_iter in xrange(0, len(barcode_haplotypes)):
        if n_seen[barcode_iter] == 0:
            barcode_not_seen.add(barcode_iter)
        elif barcode_concord[barcode_iter] / n_seen[barcode_iter] < 0.5:
            barcode_discordant.add(barcode_iter)

    logging.debug("Removed {} barcodes due to discordance".format(
            len(barcode_discordant)))
    logging.debug("Removed {} barcodes due to not seen".format(
            len(barcode_not_seen)))
    skip_barcode_indices = barcode_not_seen | barcode_discordant

    return determine_mosaicism(barcode_haplotypes, skip_barcode_indices,
                               barcode_index, variant_id, variant_barcodes)


def construct_germline_barcode_matrix(variant_barcodes,
                                      dtype=np.dtype("int32")):
    # all barcodes that have to do with any of the variants in the region?
    '''
    Construct a sparse matrix of the variants \
    supported by the barcode information.
    '''
    # Find the size of the matrix   #
    # Index the matrix rows/columns #
    variant_index = {}
    barcode_index = {}
    variant_counter = 0
    barcode_counter = 0
    for variant, properties in variant_barcodes.items():
        if not properties["is_germline"]:  # Skip non-het germline vars
            continue
        variant_index[variant] = variant_counter
        variant_counter += 1
        for allele, barcode_list in properties.items():
            if not (type(allele) is int and type(barcode_list) is list):
                continue
            for barcode in barcode_list:
                if barcode not in barcode_index:
                    barcode_index[barcode] = barcode_counter
                    barcode_counter += 1

    variant_barcode_matrix = np.zeros( #0's for no information
            (variant_counter, barcode_counter), dtype=np.dtype("int32"))

    # Populate the matrix #
    for variant, properties in variant_barcodes.items():
        if not properties["is_germline"]:
            continue
        for allele, barcode_list in properties.items():
            if not (type(allele) is int and type(barcode_list) is list):
                continue
            for barcode in barcode_list:
                #1's and 2's if provides information on allele 0 or 1
                variant_barcode_matrix[variant_index[variant],
                                       barcode_index[barcode]] = allele + 1
    return (variant_barcode_matrix, variant_index, barcode_index)


def get_haplotypes_DP(variant_matrix,ploidy=2):
    '''
    Determine the haplotypes of the barcodes.
    barcode_haplotypes: Haplotype (1 or 0) of each barcode
    barcode_concord: Fraction of variants that are concordant on each barcode
    variant_haplotypes: Allele (1 or 0) of each variant on haplotype 0 (this "haplotype 0" is numbered wrt barcode_haplotypes)
    variant_concord: fraction of barcodes concordant with this variant
    n_seen: nubmer of variants on each barcode
    '''
    if ploidy == 2:
        return get_haplotypes_diploid(variant_matrix)
    else:
        return get_haplotypes_polyploid(variant_matrix,ploidy)

def get_haplotypes_polyploid(variant_matrix,ploidy):
    #rows = barcodes; columns = variants
    variant_matrix = variant_matrix.T
    n_variants = variant_matrix.shape[1]
    n_barcodes = variant_matrix.shape[0]
    
    #initialization values for first <ploidy> cells in traceback
    #list of lists: the ith list is the case for assigning to the ith haplotype
    H = []
    C = []
    T = []
    for i in xrange(ploidy):
        H.append([1] * n_variants) #1's or 2's for alleles
        C.append([0.0] * n_variants) # concordant barcodes
        T.append([0.0] * n_variants) # nonzero barcodes 
    
    traceback = np.zeros(n_barcodes,ploidy) # index of cell that the up-arrow points to
    n_seen = [0] * n_barcodes
    haplotypes = [0] * n_barcodes
    barcode_concord = [0] * n_barcodes
    
    for i in xrange(n_barcodes):
        bc = variant_matrix[i]
        
        #temporary copies of data structures for all haplotypes
        #list (each cell in this row) of lists (possible arrow directions) 
        #of lists (assigning to each hap)
        H_tmp = [] 
        C_tmp = [] 
        T_tmp = []
        for p in xrange(ploidy): 
            H_tmp.append([vec[:] for vec in H])
            C_tmp.append([vec[:] for vec in C])
            T_tmp.append([vec[:] for vec in T])
            
        for j in xrange(n_variants): #update concordance
            if bc[j] != 0:
                n_seen[i] += 1
                for p in xrange(ploidy): #cells in row
                    for vec in T_tmp[p]: #possible arrow directions
                        vec[j] += 1 #this variant
                for p in xrange(ploidy):
                    for vec_idx, vec in H_tmp[p]:
                        if vec[j] == bc[j]:
                            C_tmp[p][vec_idx][j] += 1
                
            for p in xrange(ploidy):
                for vec_idx, vec in enumerate(H_tmp[p]):
                # Is this still the right logic?
                    if C_tmp[p][vec_idx][j]/T_tmp[p][vec_idx][j] < 0.5: #change haplotype?
                        C_tmp[p][vec_idx][j] = C_tmp[p][vec_idx][j] - T_tmp[p][vec_idx][j]
                        vec[vec_idx][j] = (1 if vec[vec_idx][j] == 2 else 2)  

        for p in xrange(ploidy): #what is the best option for each cell up arrow?
            best_idx = 0
            best_value = 0.0
            for p_possible in xrange(ploidy): #possible haplotypes for up arrow to point to
                if sum([C_tmp[p][p_possible][i]/T_tmp[p][p_possible][i] for i in xrange(n_variants)]) > best_value: 
                best_idx = p_possible
            H[p],C[p],T[p] = H_tmp[best_idx],C_tmp[best_idx],T_tmp[best_idx]
            traceback[i][p] = best_idx
     
    #what is the best in the last row?
    for p in xrange(ploidy):
        best_idx = 0
        best_value = 0.0
        if sum([C[p][i]/T[p][i] for i in xrange(n_variants)]) > best_value: 
        best_idx = p_possible
    variant_haplotypes = H[best_idx]
    variant_concord = [C[best_idx][i]/T[best_idx][i] for i in xrange(n_variants)]
    barcode_haplotypes[-1] = best_idx
    
    for i in xrange(n_barcodes): 
        bc = variant_matrix[i]
        barcode_concord[i] = 1.0*sum([1 if bc[j] == H[j] for j in xrange(n_variants)])/n_seen[i]
        
    current_col = haplotypes[-1]
    for i in xrange(2,n_barcodes+1):
        barcode_haplotypes[-i] = traceback[-i][current_col]
        current_col = traceback[-i][current_col]
     
    return (barcode_haplotypes, barcode_concord, variant_haplotypes, variant_concord, n_seen)
        
def get_haplotypes_diploid(variant_matrix):    

    #rows = barcodes; columns = variants
    variant_matrix = variant_matrix.T
    n_variants = variant_matrix.shape[1]
    n_barcodes = variant_matrix.shape[0]
    
    #initialization values for first 2 cells
    H_yes = [1] * n_variants #1's or 2's for alleles
    C_yes = [0.0] * n_variants  # concordant barcodes
    T_yes = [0.0] * n_variants # nonzero barcodes
    H_no = [1] * n_variants 
    C_no = [0.0] * n_variants  
    T_no = [0.0] * n_variants 
    
    traceback = np.zeros(n_barcodes,2) # 0 = UP; 1 = DIAG
    n_seen = [0] * n_barcodes
    haplotypes = [0] * n_barcodes
    barcode_concord = [0] * n_barcodes
    
    for i in xrange(n_barcodes):
        bc = variant_matrix[i]
        
        #case 1: this barcode isn't flipped
        H_yes_tmp1 = H_yes[:]
        C_yes_tmp1 = C_yes[:]
        T_yes_tmp1 = T_yes[:]
        H_no_tmp1 = H_no[:]
        C_no_tmp1 = C_no[:]
        T_no_tmp1 = T_no[:]
        
        #case 2: this barcode is flipped
        H_yes_tmp2 = H_yes[:]
        C_yes_tmp2 = C_yes[:]
        T_yes_tmp2 = T_yes[:]
        H_no_tmp2 = H_no[:]
        C_no_tmp2 = C_no[:]
        T_no_tmp2 = T_no[:]
        
        for j in xrange(n_variants): #update concordance
            if bc[j] != 0:
                n_seen[i] += 1
                T_yes_tmp1 += 1
                T_no_tmp1 += 1
                T_yes_tmp2 += 1
                T_no_tmp2 += 1
                if bc[j] == H_yes_tmp1[j]: #both 1's or both 2's
                    C_yes_tmp1 += 1
                if bc[j] == H_no_tmp1[j]:
                    C_no_tmp1 += 1
                if bc[j] != H_yes_tmp2[j]: #opposite
                    C_yes_tmp2 += 1
                if bc[j] != H_no_tmp2[j]:
                    C_no_tmp2 += 1
            if C_yes_tmp1[j]/T_yes_tmp1[j] < 0.5: #change haplotype?
                C_yes_tmp1[j] = T_yes_tmp1[j] - C_yes_tmp1[j]
                H_yes_tmp1[j] = (1 if H_yes_tmp1[j] == 2 else 2)
            if C_yes_tmp2[j]/T_yes_tmp2[j] < 0.5:
                C_yes_tmp2[j] = T_yes_tmp2[j] - C_yes_tmp2[j]
                H_yes_tmp2[j] = (1 if H_yes_tmp2[j] == 2 else 2) 
                
        #case 1: not flipped
        if sum([C_yes_tmp1[i]/T_yes_tmp1[i] for i in xrange(n_variants)]) > sum([C_no_tmp1[i]/T_no_tmp1[i] for i in xrange(n_variants)]): 
            H_no,C_no,T_no = H_yes_tmp1,C_yes_tmp1,T_yes_tmp1
            traceback[i][1] = 1 #DIAG
        if sum([C_yes_tmp1[i]/T_yes_tmp1[i] for i in xrange(n_variants)]) <= sum([C_no_tmp1[i]/T_no_tmp1[i] for i in xrange(n_variants)]): 
            H_no,C_no,T_no = H_no_tmp1,C_yes_tmp1,T_yes_tmp1
            traceback[i][1] = 0 #UP
        #case 2: flipped       
        if sum([C_yes_tmp2[i]/T_yes_tmp2[i] for i in xrange(n_variants)]) > sum([C_no_tmp2[i]/T_no_tmp2[i] for i in xrange(n_variants)]): 
            H_yes,C_yes,T_yes = H_yes_tmp2,C_yes_tmp2,T_yes_tmp2
            traceback[i][2] = 0 #UP
        if sum([C_yes_tmp2[i]/T_yes_tmp2[i] for i in xrange(n_variants)]) <= sum([C_no_tmp2[i]/T_no_tmp2[i] for i in xrange(n_variants)]): 
            H_yes,C_yes,T_yes = H_no_tmp2,C_yes_tmp2,T_yes_tmp2
            traceback[i][2] = 1 #DIAG
    
    #What to choose in last row
    if sum([C_yes[i]/T_yes[i] for i in xrange(n_variants)]) > sum([C_no[i]/T_no[i] for i in xrange(n_variants)]):
        variant_haplotypes = H_yes
        variant_concord = [C_yes[i]/T_yes[i] for i in xrange(n_variants)]
        barcode_haplotypes[-1] = 1
    else:
        variant_haplotypes = H_no
        variant_concord = [C_no[i]/T_no[i] for i in xrange(n_variants)]
        barcode_haplotypes[-1] = 0      
        
    for i in xrange(n_barcodes): 
        bc = variant_matrix[i]
        barcode_concord[i] = 1.0*sum([1 if bc[j] == H[j] for j in xrange(n_variants)])/n_seen[i]
    
    current_col = haplotypes[-1]
    for i in xrange(2,n_barcodes+1):
        if current_col == 0: #not-flipped column
            if traceback[-i][current_col] == 0: #stay in not-flipped state
                barcode_haplotypes[-i] = 0
            else: #went diagonally to flipped state
                barcode_haplotypes[-i] = 1
                current_col = 1
        else: #flipped column
            if traceback[-i][current_col] == 0: #stay in flipped state
                barcode_haplotypes[-i] = 1
            else: #went diagonally to not-flipped state
                barcode_haplotypes[-i] = 0
                current_col = 0
        
    return (barcode_haplotypes, barcode_concord, variant_haplotypes, variant_concord, n_seen)


def determine_mosaicism(haplotypes, skip_barcode_indices,
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
