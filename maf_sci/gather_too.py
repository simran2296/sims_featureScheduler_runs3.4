import pandas as pd
import os
import glob


if __name__ == "__main__":

    files = glob.glob("h5s/*.h5")

    table_list = [pd.read_hdf(filename, 'too') for filename in files]
    final_table = pd.concat(table_list)
    final_table.to_hdf('summary_too.h5', 'too')
