import os, sys, subprocess, time, argparse

parser = argparse.ArgumentParser("Runs the parametrized ABX example")
parser.add_argument("--maxN", default=5, type=int,
                    help="Maximal number of initial tokens")
parser.add_argument("--pn2MAXHDA", default="./build/pn2MAXHDA",
                    type=str, help="path to pn2MAXHDA executable.")
parser.add_argument("--pnmlPath", default="./pnml_files/basic/",
                     type=str, help="path to pnml files. "
                                    "Needs read and write.")
parser.add_argument("--outPath", default="./tmpData",
                    type=str, help="Where to store the csv's")

if __name__ == "__main__":

    args = parser.parse_args()

    allConf = ["standard",
               "standardmaxhda_dfs", "standardmaxhda_bfs",
               "standardmaxhda2_dfs", "standardmaxhda2_bfs"]

    allConfOut = [os.path.join(args.outPath, f"out_{ac}.csv")
                  for ac in allConf]

    os.makedirs(args.outPath, exist_ok=True)

    with open(os.path.join(args.outPath, "instance.txt"), "w") as fInst:
        fInst.write("inst\n")

        for ac, acf in zip(allConf, allConfOut):
            subprocess.call(f"""rm -f {acf}""", shell=True)

        for ll in range(1,args.maxN+1):
            for rr in range(ll, args.maxN+1):
                print("Config is ", ll, "and", rr)
                fInst.write(f"ab_{ll}_{rr}\n")

                opt = "csv_time_header" if (ll == 1) and (rr == 1) else "csv_time"

                subprocess.call(f"""cat {os.path.join(args.pnmlPath, 'abx_par.pnml')} | sed 's|<text>INITLEFT</text>|<text>{ll}</text>|g' | sed 's|<text>INITRIGHT</text>|<text>{rr}</text>|g' > {os.path.join(args.pnmlPath, 'abx_ll_rr.pnml')}""", shell=True)

                for ac, acf in zip(allConf, allConfOut):
                    print(ac, time.time())
                    subprocess.call(f"""{args.pn2MAXHDA} {ac} {os.path.join(args.pnmlPath, 'abx_ll_rr.pnml')} {opt} >> {acf}""", shell=True)
