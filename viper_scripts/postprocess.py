# Maps Ensembl IDs in a VIPER activities matrix to Hugo
# writes out translated activity matrix

import pandas as pd
import numpy as np
import argparse

def map_IDs(in_file):
    mapfile = '/home/exacloud/lustre1/CompBio/manningh/PyCharmProjects/smmartVIPER/data/target_id.txt'
    mapdf = pd.read_csv(mapfile, header=0, usecols=['V2','V5','V6'], sep=" ")
    mapdf = mapdf.set_index('V2')

    ensembl_df = pd.read_csv(in_file, header=0, index_col=0, sep='\t')

    # todo: this won't work if any don't map uniquely
    #id_map = {eid: np.unique(mapdf.loc[eid]['V6']) for eid in ensembl_df.index}

    new_index = [np.unique(mapdf.loc[eid]['V6']) for eid in ensembl_df.index]

    for i in range(len(new_index)):
        if len(new_index[i]) == 1:
            new_index[i] = new_index[i][0]
        elif len(new_index[i]) == 0:
            new_index[i] = 'UNMAPPED'
        elif len(new_index[i]) >> 1:
            new_index[i] = ", ".join(new_index[i])

    hugo_df = ensembl_df.copy()
    hugo_df.index = new_index

    return hugo_df


def main():
    print("Postprocessing...")
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory",
                        type=str,
                        help='Full path to directory containing activity matrix with Ensembl IDs in 1st col\n'
                             'Will look for file named: viper_activities_ensembl.tsv'
                        )

    args = parser.parse_args()

    ensembl_file = args.directory + 'viper_activities_ensembl.tsv'

    hugo_file = ensembl_file.replace("ensembl.tsv", "symbol.tsv")
    hugo_df = map_IDs(ensembl_file)

    hugo_df.to_csv(hugo_file, sep='\t')


if __name__ == '__main__':
    main()





