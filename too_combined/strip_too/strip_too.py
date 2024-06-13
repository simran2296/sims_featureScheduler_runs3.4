import glob
import pandas as pd
import sqlite3
import os


if __name__ == "__main__":

    files = glob.glob('../*10yrs.db')

    qu = "select * from observations where note not like 'ToO%';"
    for filename in files:
        con = sqlite3.connect(filename)
        df = pd.read_sql(qu, con)
        con.close()
        new_file = os.path.basename(filename).replace('10yrs.db', 'stripped_10yrs.db')
        con = sqlite3.connect(new_file)
        df.to_sql("observations", con)
        con.close()
