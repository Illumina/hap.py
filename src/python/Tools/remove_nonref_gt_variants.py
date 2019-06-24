#!/illumina/development/haplocompare/hc-virtualenv/bin/python
import sys
import re

GT_SPLITTER = re.compile(r'[\/\|]') # split a genotype field by '/' or '|'

def fast_nonref_remover(input_stream, output_stream):
    """
    Copies each line of :param:`input_stream` to :param:`output_stream`
    unless that line describes a variant where any of the sample genotypes
    include a <NON_REF> allele.
    """
    for line in input_stream:
        bad_variant = False            
        if line[0] != '#':
            split_line = line.split('\t')
            split_alt = split_line[4].split(',')
            if split_alt[-1] == '<NON_REF>':
                for sample_column in range(9, len(split_line)): # 
                    split_gt = re.split(GT_SPLITTER, split_line[sample_column].split(':')[0])
                    n_alt = len(split_alt)
                    for gt in split_gt:
                        if gt != '.' and int(gt) == n_alt:
                            bad_variant = True
                            break
                    else: # break out of both loops as soon as we find the first <NON_REF> GT
                        continue
                    break
        if not bad_variant:
            output_stream.write(line)
                
                
if __name__ == '__main__':
    fast_nonref_remover(sys.stdin, sys.stdout)