"""
Minimal preprocessing currently takes place. This step only consists of generating a mapped regulator file
for input into aracne. It also outputs the expression file in a now-tab-separated and non gzipped format for
ARACNE, although this may not be necessary.


Reads in
    (1) gzipped tatlow tcga expr file (ensgacs)
    (2) regulator list (hugo IDs)

Spits out
    (1) unzipped and properly formatted ensembl tatlow tcga expr file
    (2) regulator list (ensembl IDs)

Author: Hannah Manning
Date: April 2018

"""

import pandas as pd
import numpy as np
import argparse

def map_regulators(regulatorin, regulatorout, mapfile):
    """Takes a regulator input file with HUGO IDs (one per line) and maps them to ensembl IDs.
    Generates a new regulator file with ensembl IDs.
    """

    mapdf = pd.read_csv(mapfile, header=0, usecols=['V2', 'V5', 'V6'], sep=" ")
    mapdf = mapdf.set_index('V2')

    with open(regulatorin, "r") as fh:
        regulators = fh.read().splitlines()

    # map regulators to ensembl ids
    reg_ensgacs = []
    for r in regulators:
        ensgacs = np.unique(mapdf.index[mapdf['V6'] == r]).tolist()
        for e in ensgacs:
            reg_ensgacs.append(e)

    # eliminate those which don't appear in the expr dataframe
    # todo: this doesn't happen in the brca so i'm not worrying about it for now!
    # in future, drop those from reg_ensgacs list if they don't appear in exprdf.index

    reg_fh = open(regulatorout, "w")
    for r in reg_ensgacs:
        reg_fh.write("%s\n" % r)


def main():
    parser = argparse.ArgumentParser()

    #parser.add_argument("-ei", "--expression_in",
    #                    type=str,
    #                    help='Full path to expression file (Kallisto-processed by Tatlow et al.\n'
    #                         'i.e. path/to/TCGA_BRCA_bygene.txt.gz; (can be gzipped or not)'
    #                    )

    #parser.add_argument("-eo", "--expression_out",
    #                    type=str,
    #                    help='Full path to new expression file (filtered, mapped).'
    #                    )

    parser.add_argument("-ro", "--regulator_out",
                        type=str,
                        help='Full path to new (filtered, mapped) regulator file'
                        )

    args = parser.parse_args()

    #exprin = args.expression_in
    #exprout = args.expression_out
    regulatorout = args.regulator_out

    lustrehome = '/home/exacloud/lustre1/CompBio/manningh/'
    #regulatorin = lustrehome + 'Data/tf_act/sens_val/pancancer-tfgenes_2014_04_17_filtered.txt'
    regulatorin = lustrehome + 'Data/tf_act/sens_val/PC_regulators.txt'
    mapfile = lustrehome + 'PyCharmProjects/smmartVIPER/data/target_id.txt'

    # todo: add compression checker
    #expr = pd.read_csv(exprin, compression='gzip', sep=" ", header=0, index_col=1)
    #expr = expr.drop('Unnamed: 0', 1)

    map_regulators(regulatorin, regulatorout, mapfile)

    # todo: this step is not currently necessary, but may have to write out expr file like this if other
    # todo: preprocessing occurs
    #expr.to_csv(exprout, sep='\t')


if __name__ == '__main__':
    main()
