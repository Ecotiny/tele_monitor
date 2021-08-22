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
    print("Selecting mag between VT and BT")
    df['mag'] = df.apply(mag_selector, axis=1)
    print("Selecting only stars with mean positions")
    len_b = len(df)
    df = df[df.pflag != "X"]
    len_a = len(df)
    print(f"Removed {len_b - len_a} stars")

    print("Getting only important columns")
    df = df[["HIP", "mRAdeg", "mDEdeg", 'mag']]

    print("Converting to float positions")
    df["mRAdeg"] = df["mRAdeg"].astype(float)
    df["mDEdeg"] = df["mDEdeg"].astype(float)

    return df

def process_hip(fn, df):
    hip = pd.read_csv(
            fn,
            delimiter="|")
    hip[["RArad", "DErad"]] = hip['        RArad         DErad'].str.split(expand=True)
    hip.columns = hip.columns.str.strip() # strip whitespace from column names
    df = df[~df.HIP.isin(hip.HIP)]    
    hip['mRAdeg'] = hip['RArad'].astype(float) * 360/(2*np.pi)
    hip['mDEdeg'] = hip['DErad'].astype(float) * 360/(2*np.pi)
    hip['mag'] = hip['Hpmag'].astype(float)

    hip = hip[["HIP", "mRAdeg", "mDEdeg", 'mag']]

    df = df.append(hip)
    return df 


if __name__ == "__main__":
    overall_df = pd.DataFrame() 
    for n in range(20):
        fn = f"tyc2.dat.{n:02d}.gz"
        print(fn)
        overall_df = overall_df.append(process_file(fn))

    print('intersecting with hipparcos')
    overall_df = process_hip("hip2.dat", overall_df)
    print('Removing NaN')
    overall_df.dropna(inplace=True)
    print("Converting to SQL")

    db = sqlite3.connect("tycho2.db")
    overall_df.to_sql("catalogue", db, if_exists='replace', index=False)
    print(overall_df)
