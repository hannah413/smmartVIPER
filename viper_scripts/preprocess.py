# Incorporates a single smmart patient and a UHR into a gene expression matrix
# preprocess.py -b ${ORIG_EXPR} -n ${NEW_EXPR} -s ${SOI} -u ${UHR}

import argparse
import pandas as pd
import numpy as np


def log2_normalize(df):
    """
    This is a placeholder. Will not be used in SMMART trials as is.
    Log2-normalizes an expression df.
    Performed in same fashion as in bergamot utils.
    Not going to worry about the double-transposition as this function will be removed.

    Args:
        df
            pandas.DataFrame of floats
            shape = (n_features, n_samples)

    Returns:
        norm_df
            pandas.DataFrame of floats
            shape = (n_features, n_samples)
    """

    df = df.T
    log_add = np.nanmin(df[df > 0]) * 0.5
    norm_df = np.log2(df + log_add)
    norm_df = norm_df.T

    return norm_df


def check_name(id):
    """
    Adds 'X' in front of ID if not present.
    """

    if not id.startswith('X'):
        id = 'X' + id

    return id


def main():
    print("Preprocessing...")
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--background_expr",
                        type=str,
                        help="Full path to background expression matrices (i.e. in tatlow_matrices/)"
                        )

    parser.add_argument("-n", "--new_expr",
                        type=str,
                        help="Full path to output file of this script (processed expression)"
                        )

    parser.add_argument("-s", "--sample_of_interest",
                        type=str,
                        default=None,
                        help="ID of sample of interest in /patterja/smmart_matrices/counts_smmart_unfilt.txt\n"
                             "(ID may be chosen from ProcessedID col of md_smmart.txt in same dir)\n"
			     "i.e. X111111_NNN00000_0000_AAAA0AAAA0_RNA111111AA_DNA.11.11111"
                        )

    parser.add_argument("-u", "--universal_human_ref",
                        type=str,
                        default=None,
                        help="ID of universal human reference run with sample of interest\n"
                             "in patterja/smmart_matrices_counts_smmart_unfilt.txt\n"
			     "i.e. 111111_NN000000_0000_AAAA0AAAA0_RNA111111AA_UHR"
                        )


    args = parser.parse_args()
    orig_expr_fl = args.background_expr
    new_expr_fl = args.new_expr

    soi = args.sample_of_interest
    uhr = args.universal_human_ref

    soi = check_name(soi)
    uhr = check_name(uhr)

    orig_expr = pd.read_csv(orig_expr_fl, compression='gzip', sep=" ", header=0, index_col=1)
    orig_expr = orig_expr.drop('Unnamed: 0', 1)

    new_expr = orig_expr.copy()

    # add columns for the soi and uhr into the expression dataset
    all_smmart_expr = pd.read_csv("/home/exacloud/lustre1/CompBio/users/patterja/smmart_matrices/counts_smmart_unfilt.txt",
                                  sep=" ", header=0, index_col=0)

    # TODO: add appropriate normalization functions!
    if soi:
        new_expr[soi] = all_smmart_expr.loc[:, soi]
    if uhr:
        new_expr[uhr] = all_smmart_expr.loc[:, uhr]

    new_expr = log2_normalize(new_expr)

    # write out the preprocessed expression file
    new_expr.to_csv(new_expr_fl, sep='\t')


if __name__ == '__main__':
    main()
