import pandas as pd
import argparse
import os

def getoptions():
    parser = argparse.ArgumentParser(description='Extract counts information from lima outputs')
    parser.add_argument("-d", "--dir", dest="dir", required=True, help="Path to the directory where the input files are in")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Path to the directory where the output files will be in")
    parser.add_argument("-p", "--preffix", dest="preffix", required=True, help="A preffix for all the output files")
    args = parser.parse_args()
    return(args)


# extract information from lima standard output files *.lima.counts and summarize the counts information
def summary_from_files(path, outpath, preffix):
    path = os.path.abspath(path)
    outpath = os.path.abspath(outpath)
    files = os.listdir(path)
    files = [i for i in files if (i!="log" and i!="summary")] 
    df = pd.DataFrame()
    for i in files:
        fpath = path + "/" + i + "/" + i + "_limademuxed.lima.counts"
        ff = pd.read_table(fpath)
        df0 = {'fileName' : [i] * ff.shape[0], 'sampleID' : [q.split("_3p")[0] for q in ff['IdxCombinedNamed'].tolist()], 'lima-counts' : ff['Counts'].tolist()}
        df0 = pd.DataFrame(df0)
        df = pd.concat([df, df0], axis=0)

    df.to_csv(outpath + '/' + preffix + "_lima_summary.csv", index=False)
    df_by_sample = df.drop('fileName', axis=1).groupby("sampleID").sum()
    df_by_sample.to_csv(outpath + '/' + preffix + "_lima_sampleSum.csv")
        

def main():
    args = getoptions()
    summary_from_files(args.dir, args.out, args.preffix)
    return

if __name__=="__main__":
    main()
