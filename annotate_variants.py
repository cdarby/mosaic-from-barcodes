#!/usr/bin/env python3

import argparse
import sys
import logging

# Handle b37 and hg19 chromosomes #
chroms = dict([[str(x), x] for x in list(range(1, 23))])
chroms.update({'X': 23, 'Y': 24, "MT": 25})
chroms.update(dict([["chr" + str(x), x] for x in list(range(1, 23))]))
chroms.update({"chrX": 23, "chrY": 24, "chrM": 25})
# Add a last chromosome #
chroms["END"] = 9999

log_format = "%(filename)s::%(funcName)s [%(levelname)s] %(message)s"


def process_args():
    parser = argparse.ArgumentParser(
        description="Annotate variant present in other VCF files")
    parser.add_argument("--infile",
                        help="The input VCF file to annotate [stdin]")
    parser.add_argument("--outfile", help="The output VCF file [stdout]")
    parser.add_argument("--ref_var_file", nargs='+',
                        help="A VCF file containing external variants")
    parser.add_argument("--ref_var_name", nargs='+',
                        help="The INFO field tag to use for "
                        "variants in the external file")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def find_first_lines(*vcf_lines):
    '''
    Given a list of VCF lines (as lists), \
    return the indexs of the VCF with the earliest variant.
    '''
    first_chrom, first_pos = ("END", 99999999)
    first_idxs = []
    for vcf_idx, line in enumerate(vcf_lines):
        if (chroms[line[0]] < chroms[first_chrom] or
            (chroms[line[0]] == chroms[first_chrom] and
             int(line[1]) < int(first_pos))):
            # The variant is the first seen #
            first_idxs = [vcf_idx]
            first_chrom, first_pos = line[:2]
        elif (chroms[line[0]] == chroms[first_chrom] and
              line[1] == first_pos): #in case of repeat
            first_idxs.append(vcf_idx)
    return first_idxs


def read_next_lines(ref_fhs, ref_lines, ref_line_idx):
    '''
    Read the next VCF lines from file handles with the given indexes.
    '''
    for idx in ref_line_idx:
        if idx < 0:  # corresponds to the input VCF: its next line 
            continue # gets read elsewhere
        ref_lines[idx] = ref_fhs[idx].readline().rstrip().split('\t')
        if not ref_lines[idx] or not ref_lines[idx][0]:
            ref_lines[idx] = ["END", "9999999"]


def main(args):
    if not args:
        args = process_args()

    log_level = logging.WARNING
    if args.debug:
        log_level = logging.DEBUG
    logging.basicConfig(format=log_format, level=log_level)

    if not args.ref_var_file or not args.ref_var_name:
        sys.exit(1)
    if len(args.ref_var_file) != len(args.ref_var_name):
        sys.exit("Please supply a name for each input reference VCF")

    if args.infile:
        args.in_fh = open(args.infile)
    else:
        args.in_fh = sys.stdin

    if args.outfile:
        args.out_fh = open(args.outfile)
    else:
        args.out_fh = sys.stdout

    logging.debug("Supplied ref files: {}".format(
            str(args.ref_var_file)))
    ref_fhs = []
    for ref_file in args.ref_var_file:
        ref_fhs.append(open(ref_file))
    logging.debug("Ref file handles: {}".format(
            str(ref_fhs)))

    # Get to the first data line for all VCFs #
    # infile
    vcf_line = args.in_fh.readline().rstrip()
    while vcf_line[0] == '#':
        print(vcf_line, file=args.out_fh)
        vcf_line = args.in_fh.readline().rstrip()
    vcf_line = vcf_line.split('\t')

    # other files
    ref_lines = [''] * len(ref_fhs)
    for ref_idx, ref_fh in enumerate(ref_fhs):
        ref_lines[ref_idx] = ref_fh.readline().rstrip()
        while ref_lines[ref_idx][0] == '#':
            ref_lines[ref_idx] = ref_fh.readline().rstrip()
        ref_lines[ref_idx] = ref_lines[ref_idx].split('\t')
    logging.debug("Ref lines: {}".format(
            str(ref_lines)))

    # Read through VCFs and add annotations #
    # If a variant in infile is found in ref_var_file(s), annotate it accordingly
    while vcf_line:
        logging.debug("At {}".format(
                ','.join([str(x[0]) + ':' + str(x[1])
                          for x in [vcf_line] + ref_lines])))
        first_line_idxs = find_first_lines(vcf_line, *ref_lines)
        logging.debug("First line indexs: {}".format(
                str(first_line_idxs)))
        if 0 in first_line_idxs: # 0 indicates infile
            logging.debug("Reading new vcf line")
            for variant_idx in first_line_idxs[1:]: #duplicates
                vcf_line[7] += ';' + args.ref_var_name[variant_idx - 1]
            #print the new entry
            print('\t'.join(vcf_line), file=args.out_fh)
            #read next line in infile
            vcf_line = args.in_fh.readline().rstrip().split('\t')
            if not vcf_line or not vcf_line[0]:
                vcf_line = None
        read_next_lines(ref_fhs, ref_lines, [x - 1 for x in first_line_idxs]) 
        #get next line of VCFs that were used to annotate this line in infile

    if args.infile:
        args.in_fh.close()
    if args.outfile:
        args.out_fh.close()
    for ref_fh in ref_fhs:
        ref_fh.close()

if __name__ == "__main__":
    main(None)
