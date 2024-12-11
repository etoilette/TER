import os, sys, subprocess, argparse, tempfile

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import pandas

def myplot(csv: str):
    x = pandas.read_csv(csv)

    N = x.shape[0] - 1

    colors = [cm.jet(c) for c in np.linspace(0, 1, N)]
    ff,aa = plt.subplots(1,1)
    for idx, l in enumerate(x.iterrows()):
        d = np.array(eval(l[1]["L(dim;ncell)"].replace(";", ","))).T
        red = np.array(eval(l[1]["MaxCells"].replace(";", ","))).T
        name = str(idx)
        try:
            name = l[1]["name"]
        except:
            pass
        aa.plot(d[0,:], d[1,:], ".-", markersize=10, color=colors[idx%len(colors)], label=name + " All Cells")
        aa.plot(red[0,:], red[1,:], "*", markersize=10, color=colors[idx%len(colors)], label=name + " Max Cells")

    aa.legend()
    aa.grid("on")
    aa.set_yscale("log")
    aa.set_xlabel("conclist dim")
    aa.set_ylabel("# of cells")
    plt.show()

def run_statistics(args, fold: str, tnames, tcsv):
    """
    Runs statistics on all pnml files found in fold and the subfolders
    """
    for f in os.listdir(fold):
        if os.path.isfile(os.path.join(fold, f)):
            if f.endswith(".pnml"):
                n = os.path.splitext(f)[0]
                print(f"Treating {n}")
                tnames.write(f"{n}\n")
                subprocess.check_call(f""" {args.pn2HDA} "standard" {os.path.join(fold, f)} "csv_header" >> {tcsv.name}""",
                                      shell=True)
            else:
                continue
        else:
            run_statistics(args, os.path.join(fold, f), tnames, tcsv)

ap = argparse.ArgumentParser(prog="show_pHDA", usage="none", description="Constructs a (partial) HDA from pnml files "
                                                                         "and outputs its statistics.")
ap.add_argument("--pn2HDA", default="./build/pn2HDA", type=str, help="path to pn2HDA executable.")
ap.add_argument("--csv", type=str, default="./examples/res.csv",
                help="Name of the csv to generate or use.")
ap.add_argument("--runstats", action="store_true",
                help="If set, computes the examples first, otherwise reuses the existing csv file.")
ap.add_argument("--instances", type=str, default="pnml_files/mcc/")

if __name__ == "__main__":
    args = ap.parse_args()

    if args.runstats:
        with tempfile.NamedTemporaryFile("w") as tnames, tempfile.NamedTemporaryFile("w") as tcsv:
            tnames.write("name\n")
            run_statistics(args, args.instances, tnames, tcsv)
            os.makedirs(os.path.dirname(args.csv), exist_ok=True)
            tnames.flush()
            tcsv.flush()
            with tempfile.NamedTemporaryFile("w") as tcsv2:
                subprocess.check_call(f"head -n1 {tcsv.name} > {tcsv2.name}", shell=True)
                subprocess.check_call(f"""awk 'NR % 2 == 0' {tcsv.name} >> {tcsv2.name}""", shell=True)
                subprocess.call(f"""paste -d "," {tnames.name} {tcsv2.name} > {args.csv}""", shell=True)

    myplot(args.csv)
