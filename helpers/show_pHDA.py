import sys, os, subprocess, argparse


ap = argparse.ArgumentParser(prog="show_pHDA", usage="none", description="Constructs ST-automata from pnml files "
                                                                         "and creates the corresponding pdfs.")
ap.add_argument("pnml_file", type=str, help="pnml file or folder containing pnml files to process."
                                                         "If it is a folder, then all files will be converted.")
ap.add_argument("--out_dir", type=str, default="./examples/",
                help="output directory for dot and pdf files.")
ap.add_argument("--pn2HDA", default="./build/pn2HDA", type=str, help="path to pn2HDA executable.")

def gen1(args, pnml: str, out_dir: str, show:bool):
    n = os.path.splitext(os.path.basename(pnml))[0]
    subprocess.check_call(f"""{args.pn2HDA} "wpni" "{pnml}" "dot" > {os.path.join(out_dir, n + ".dot")}""", shell=True)
    subprocess.check_call(["dot", "-Tpdf", "-o", os.path.join(out_dir, n + ".pdf"), os.path.join(out_dir, n + ".dot")])

    if show:
        subprocess.call(f"open {os.path.join(out_dir, n + '.pdf')}", shell=True)


if __name__ == '__main__':
    args = ap.parse_args()

    subprocess.check_output(["mkdir", "-p", args.out_dir])

    if os.path.isfile(args.pnml_file):
        gen1(args, args.pnml_file, args.out_dir, True)
    else:
        for f in os.listdir(args.pnml_file):
            if not f.endswith(".pnml"):
                continue
            gen1(args, os.path.join(args.pnml_file, f), args.out_dir, False)

    sys.exit(0)

