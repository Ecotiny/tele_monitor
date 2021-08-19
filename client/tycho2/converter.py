import sqlite3
import os
import numpy as np
import pandas as pd

def mag_selector(row):
    if str(row['VT']).strip() != "":
        return float(row['VT'])
    return float(row['BT'])

def process_file(fn):
    print("Importing")
    df = pd.read_csv(
            fn,
            compression="gzip",
            delimiter="|",
            names=["CATNUMS", "pflag", 
                   "mRAdeg", "mDEdeg", "pmRA", "pmDE",
                   "e_mRA", "e_mDE", "e_pmRA", "e_pmDE",
                   "mepRA", "mepDE", "Num", 
                   "g_mRA", "g_mDE", "g_pmRA", "g_pmDE",
                   "BT", "e_BT", "VT", "e_VT", "prox",
                   "TYC", "HIP", 
                   "RAdeg", "DEdeg", "epRA", "epDE",
                   "e_RA", "e_DE", "poslfg", "corr"])
    print("Generating catalogue columns")
    df[["TYC1", "TYC2", "TYC3"]] = df.CATNUMS.str.split(" ", expand=True)
    print("Dropping catnums")
    df = df.drop("CATNUMS", axis="columns")
    print("Selecting only stars with mean positions")
    df = df[df.pflag != "X"]
    print("Selecting mag between VT and BT")
    df['mag'] = df.apply(mag_selector, axis=1)

    print("Getting only important columns")
    df = df[["TYC2", "mRAdeg", "mDEdeg", 'mag']]

    print("Converting to float positions")
    df["mRAdeg"] = df["mRAdeg"].astype(float)
    df["mDEdeg"] = df["mDEdeg"].astype(float)

    return df


if __name__ == "__main__":
    overall_df = pd.DataFrame() 
    for n in range(20):
        fn = f"tyc2.dat.{n:02d}.gz"
        print(fn)
        overall_df = overall_df.append(process_file(fn))

    print("Converting to SQL")

    db = sqlite3.connect("tycho2.db")
    overall_df.to_sql("catalogue", db, if_exists='replace', index=False)
    print(overall_df)
