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
    #(barcode_haplotypes, barcode_concord, n_seen) = get_haplotypes_DP(variant_matrix)
    
    skip_barcode_indices = set()
    barcode_not_seen = set()
    barcode_discordant = set()
    for barcode_iter in range(0, len(barcode_haplotypes)):
        if n_seen[barcode_iter] == 0:
            barcode_not_seen.add(barcode_iter)
        elif barcode_concord[barcode_iter] / n_seen[barcode_iter] < 0.5:
            barcode_discordant.add(barcode_iter)

    logging.debug("Removed {} barcodes due to discordance".format(
            len(barcode_discordant)))
    logging.debug("Removed {} barcodes due to not seen".format(
            len(barcode_not_seen)))
    skip_barcode_indices = barcode_not_seen | barcode_discordant

    (depth, prob_mosaic, prob_false_positive) = determine_mosaicism(barcode_haplotypes, skip_barcode_indices,
                               barcode_index, variant_id, variant_barcodes)
    
    #if prob_mosaic == 0.98:
        #print(str(barcode_concord))
        #print(str(variant_haplotypes))
    '''
    print("Removed {} barcodes due to discordance".format(
        len(barcode_discordant)))
    print("Removed {} barcodes due to not seen".format(
        len(barcode_not_seen)))
    '''
    return (depth, prob_mosaic, prob_false_positive)

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

   
def get_haplotypes_diploid(variant_matrix):    
    '''
    Arguments    
    variant_matrix: values are 0/1/2
    rows = barcodes; columns = variants
    
    Return Values
    barcode_haplotypes: values are 0/1 repr haplotypes
    variant_haplotypes: values are 0/1 repr alleles on 
        haplotype 0 from barcode_haplotypes
    barcode_concord: values are counts of concordant 
        Noticed that above, we look at concordance/num_seen < 0.5 to discard
        so it doesn't seem right to do the add/subtract method
    variant_concord: fraction of barcodes that are concordant
    n_seen: number of nonzero values from each barcode
    '''
   
    variant_matrix = variant_matrix.T
    n_barcodes = variant_matrix.shape[0]
    n_variants = variant_matrix.shape[1]
    #print(str(n_variants) + "," + str(n_barcodes))
    if n_variants < 1 or n_barcodes < 1:
        return None
    
    #initialization values for first 2 cells (in each column of DP matrix)
    H_yes = [1] * n_variants #represents 1's or 2's for alleles
    C_yes = [0.0] * n_variants  # concordant barcodes
    T_yes = [0.0] * n_variants # nonzero barcodes
    H_no = [1] * n_variants 
    C_no = [0.0] * n_variants  
    T_no = [0.0] * n_variants 
    
    traceback = np.zeros((n_barcodes,2),dtype=np.dtype("int32")) # 0 = UP; 1 = DIAG
    n_seen = [0] * n_barcodes
    barcode_haplotypes = [0] * n_barcodes
    barcode_concord = [0] * n_barcodes
    variant_haplotypes = [0] * n_variants

    curr_bc = 0
    stop = 0
    nonzero = variant_matrix.nonzero()
    while stop < len(nonzero[0]):
        #iterate through barcodes
        
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
        
        #net changes to concordance vector for each of the 4 cases considered
        sum1yes = 0.0
        sum1no = 0.0
        sum2yes = 0.0
        sum2no = 0.0
        
        while stop < len(nonzero[0]) and nonzero[0][stop] == curr_bc:
            #iterate through variants
            j = nonzero[1][stop] #variant index
            n_seen[curr_bc] += 1
            a = variant_matrix[curr_bc,j] #1 or 2 (allele)
            
            #increment count
            T_yes_tmp1[j] += 1
            T_no_tmp1[j] += 1
            T_yes_tmp2[j] += 1
            T_no_tmp2[j] += 1
            
            #compare alleles
            #want equal-to if the barcode is not flipped
            if a == H_yes_tmp1[j]: #both 1's or both 2's
                C_yes_tmp1[j] += 1
            if a == H_no_tmp1[j]:
                C_no_tmp1[j] += 1
            #want unequal-to if the barcode is flipped
            if a != H_yes_tmp2[j]:
                C_yes_tmp2[j] += 1
            if a != H_no_tmp2[j]:
                C_no_tmp2[j] += 1
            
            stop += 1    
            
            #haplotype for this variant changes if #concordant/#total is less than 0.5
            if T_yes_tmp1[j] != 0 and C_yes_tmp1[j]/T_yes_tmp1[j] < 0.5:
                C_yes_tmp1[j] = T_yes_tmp1[j] - C_yes_tmp1[j]
                H_yes_tmp1[j] = (1 if H_yes_tmp1[j] == 2 else 2)
            elif T_yes_tmp1[j] != 0 and C_yes_tmp1[j]/T_yes_tmp1[j] == 0.5:
                H_yes_tmp1[j] = 1
                #print("tie")
            sum1yes += (C_yes_tmp1[j]/T_yes_tmp1[j] if T_yes_tmp1[j] != 0 else 0.0)
            
            if T_no_tmp1[j] != 0 and C_no_tmp1[j]/T_no_tmp1[j] < 0.5:
                C_no_tmp1[j] = T_no_tmp1[j] - C_no_tmp1[j]
                H_no_tmp1[j] = (1 if H_no_tmp1[j] == 2 else 2)
            elif T_no_tmp1[j] != 0 and C_no_tmp1[j]/T_no_tmp1[j] == 0.5:
                H_no_tmp1[j] = 1
                #print("tie")
            sum1no += (C_no_tmp1[j]/T_no_tmp1[j] if T_no_tmp1[j] != 0 else 0.0)
            
            if T_yes_tmp2[j] != 0 and C_yes_tmp2[j]/T_yes_tmp2[j] < 0.5:
                C_yes_tmp2[j] = T_yes_tmp2[j] - C_yes_tmp2[j]
                H_yes_tmp2[j] = (1 if H_yes_tmp2[j] == 2 else 2) 
            elif T_yes_tmp2[j] != 0 and C_yes_tmp2[j]/T_yes_tmp2[j] == 0.5:
                H_yes_tmp2[j] = 1
                #print("tie")
            sum2yes += (C_yes_tmp2[j]/T_yes_tmp2[j] if T_yes_tmp2[j] != 0 else 0.0)
            
            if T_no_tmp2[j] != 0 and C_no_tmp2[j]/T_no_tmp2[j] < 0.5:
                C_no_tmp2[j] = T_no_tmp2[j] - C_no_tmp2[j]
                H_no_tmp2[j] = (1 if H_no_tmp2[j] == 2 else 2) 
            elif T_no_tmp2[j] != 0 and C_no_tmp2[j]/T_no_tmp2[j] == 0.5:
                H_no_tmp2[j] = 1
                #print("tie")
            sum2no += (C_no_tmp2[j]/T_no_tmp2[j] if T_no_tmp2[j] != 0 else 0.0)
        
        #determine which case was better in the traceback for each cell
        #if change from yes col to no col, it's a 1; else 0
        
        #case 1: not flipped
        if sum1yes > sum1no:
            H_no,C_no,T_no = H_yes_tmp1,C_yes_tmp1,T_yes_tmp1
            traceback[curr_bc][0] = 1 #DIAG
        else: 
            H_no,C_no,T_no = H_no_tmp1,C_no_tmp1,T_no_tmp1
            traceback[curr_bc][0] = 0 #UP
        #if H_no_tmp1 != H_yes_tmp1:
            #print("traceback tie")
        #case 2: flipped       
        if sum2yes > sum2no:
            H_yes,C_yes,T_yes = H_yes_tmp2,C_yes_tmp2,T_yes_tmp2
            traceback[curr_bc][1] = 0 #UP
        else:
            H_yes,C_yes,T_yes = H_no_tmp2,C_no_tmp2,T_no_tmp2
            traceback[curr_bc][1] = 1 #DIAG
        #if sum2yes == sum2no:
            #print("traceback tie")
            
        curr_bc += 1
    
    #What to choose in last row
    yes = [C_yes[i]/T_yes[i] if T_yes[i] != 0 else 0.0 for i in range(n_variants)]
    no = [C_no[i]/T_no[i] if T_no[i] != 0 else 0.0 for i in range(n_variants)]

    if sum(yes) > sum(no):
        variant_haplotypes = [h-1 for h in H_yes] #correct to be values 0/1
        variant_concord = yes
        barcode_haplotypes[-1] = 1
    else:
        variant_haplotypes = [h-1 for h in H_no]
        variant_concord = no
        barcode_haplotypes[-1] = 0      
    if sum(yes) == sum(no):
        print("final row tie")
        
    #construct barcode haplotypes from traceback
    current_col = barcode_haplotypes[-1]
    for i in range(2,n_barcodes+1):
        if current_col == 0: #not-flipped column
            if traceback[-i+1][current_col] == 0: #stay in not-flipped state
                barcode_haplotypes[-i] = 0
            else: #went diagonally to flipped state
                barcode_haplotypes[-i] = 1
                current_col = 1
        else: #flipped column
            if traceback[-i+1][current_col] == 0: #stay in flipped state
                barcode_haplotypes[-i] = 1
            else: #went diagonally to not-flipped state
                barcode_haplotypes[-i] = 0
                current_col = 0
    
    #construct barcode concordance from the barcode and variant haplotypes
    curr_bc = 0
    stop = 0
    while stop < len(nonzero[0]):
        while stop < len(nonzero[0]) and nonzero[0][stop] == curr_bc:
            j = nonzero[1][stop] #col in var mx
            a = variant_matrix[curr_bc,j] - 1 #convert value of 1/2 to value of 0/1
            if barcode_haplotypes[curr_bc] == 0: #not-flipped state
                if a == variant_haplotypes[j]:
                    barcode_concord[curr_bc] += 1
            else:
                if a != variant_haplotypes[j]: #flipped state
                    barcode_concord[curr_bc] += 1
            stop += 1   
        curr_bc += 1

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
                print("Constructed phasing matrix {}".format(
                    str(phasing_evidence)))
                print(variant_id)
                print(str(len(skip_barcode_indices)))
            if phasing_evidence[abs(mosaic_haplotype - 1), 1] > 2:
                prob_false_positive = 0.98

    return (depth, prob_mosaic, prob_false_positive)

    '''
    Ignore for now, just sketching out an idea
    '''
def get_haplotypes_polyploid(variant_matrix,ploidy):
    #rows = barcodes; columns = variants
    variant_matrix = variant_matrix.T
    n_variants = variant_matrix.shape[1]
    n_barcodes = variant_matrix.shape[0]
   
    #print(str(n_barcodes) + "," + str(n_variants))
    #initialization values for first <ploidy> cells in traceback
    #list of lists: the ith list is the case for assigning to the ith haplotype
    H = []
    C = []
    T = []
    for i in range(ploidy):
        H.append([1] * n_variants) #1's or 2's for alleles
        C.append([0.0] * n_variants) # concordant barcodes
        T.append([0.0] * n_variants) # nonzero barcodes 
    
    traceback = np.zeros((n_barcodes,ploidy),dtype=np.dtype("int32")) # index of cell that the up-arrow points to
    n_seen = [0] * n_barcodes
    barcode_haplotypes = [0] * n_barcodes
    barcode_concord = [0] * n_barcodes
    variant_haplotypes = [0] * n_variants
    
    for i in range(n_barcodes):
        bc = variant_matrix[i]
        
        #temporary copies of data structures for all haplotypes
        #list (each cell in this row) of lists (possible arrow directions) 
        #of lists (assigning to each hap)
        H_tmp = [] 
        C_tmp = [] 
        T_tmp = []
        for p in range(ploidy): 
            H_tmp.append([vec[:] for vec in H])
            C_tmp.append([vec[:] for vec in C])
            T_tmp.append([vec[:] for vec in T])
            
        for j in range(n_variants): #update concordance
            if bc[j] != 0:
                n_seen[i] += 1
                for p in range(ploidy): #cells in row
                    for vec in T_tmp[p]: #possible arrow directions
                        vec[j] += 1 #this variant
                for p in range(ploidy):
                    for vec_idx, vec in H_tmp[p]:
                        if vec[j] == bc[j]:
                            C_tmp[p][vec_idx][j] += 1
                
            for p in range(ploidy):
                for vec_idx, vec in enumerate(H_tmp[p]):
                # Is this still the right logic?
                    if C_tmp[p][vec_idx][j]/T_tmp[p][vec_idx][j] < 0.5: #change haplotype?
                        C_tmp[p][vec_idx][j] = C_tmp[p][vec_idx][j] - T_tmp[p][vec_idx][j]
                        vec[vec_idx][j] = (1 if vec[vec_idx][j] == 2 else 2)  

        for p in range(ploidy): #what is the best option for each cell up arrow?
            best_idx = 0
            best_value = 0.0
            for p_possible in range(ploidy): #possible haplotypes for up arrow to point to
                if sum([C_tmp[p][p_possible][i]/T_tmp[p][p_possible][i] for i in range(n_variants)]) > best_value: 
                    best_idx = p_possible
            H[p],C[p],T[p] = H_tmp[best_idx],C_tmp[best_idx],T_tmp[best_idx]
            traceback[i][p] = best_idx
    return (get_haplotypes(variant_matrix)) 
    #what is the best in the last row?
    for p in range(ploidy):
        best_idx = 0
        best_value = 0.0
        if sum([C[p][i]/T[p][i] for i in range(n_variants)]) > best_value: 
            best_idx = p_possible
    variant_haplotypes = H[best_idx]
    variant_concord = [C[best_idx][i]/T[best_idx][i] for i in range(n_variants)]
    barcode_haplotypes[-1] = best_idx
    
    for i in range(n_barcodes): 
        bc = variant_matrix[i]
        barcode_concord[i] = 1.0*sum([bc[j] == H[j] for j in range(n_variants)])/n_seen[i]
        
    current_col = haplotypes[-1]
    for i in range(2,n_barcodes+1):
        barcode_haplotypes[-i] = traceback[-i][current_col]
        current_col = traceback[-i][current_col]
     
    return (barcode_haplotypes, barcode_concord, variant_haplotypes, variant_concord, n_seen)
