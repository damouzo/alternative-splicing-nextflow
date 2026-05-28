import sys
import os
import argparse
from . import config


def getSGlist(outdir, signature_name):
    fout = open(outdir + "/" + signature_name + "/" + signature_name + "_sig_list.txt", "w")
    keys = {}
    n = 0
    for line in open(outdir + "/" + signature_name + "/" + signature_name + "_high_cor_matrix.txt"):
        if n == 0:
            n += 1
            continue
        ls = line.strip().split("\t")
        keys[ls[0].split("_")[0]] = ""
    for k in keys:
        fout.write(k + "\n")
    fout.close()


def getBGlist(outdir, signature_name):
    bg_list_name = outdir + "/" + signature_name + "/" + signature_name + "_background_list.txt"
    fout = open(bg_list_name, "w")
    keys = {}
    n = 0
    for line in open(outdir + "/" + signature_name + "/" + signature_name + "_global_cor_matrix.txt"):
        if n == 0:
            n += 1
            continue
        ls = line.strip().split("\t")
        keys[ls[0]] = ""
    for k in keys:
        fout.write(k + "\n")
    fout.close()
    return bg_list_name


def main(args):
    fin_signature    = args.signatureScorebySample
    signature_name   = fin_signature.split("/")[-1].split(".scores.txt")[0]
    fin_PSI_matrix   = args.PSIbySample
    sample_order_fin = args.groupNameOrder
    GO_analysis      = args.GO
    GO_analysis_only = args.GO_only
    if GO_analysis_only:
        bg_list = args.GO_background_gene_list

    outdir = args.out_dir

    print("###Correlation analysis")
    if not GO_analysis_only:
        # Step 1: sort / reorder signature scores by group
        cmd1 = sys.executable + " " + config.MAT_REORDER + " " + fin_signature + " " + outdir + " " + sample_order_fin
        print(cmd1)
        os.system(cmd1)

        # Step 2: generate PSI matrix aligned to sample order
        cmd2 = (
            sys.executable + " " + config.MAT_GENERATE
            + " " + fin_PSI_matrix
            + " " + outdir + "/" + signature_name + "/" + signature_name + ".sorted.txt"
        )
        print(cmd2)
        os.system(cmd2)

        # Step 3: Pearson correlation with permutation test (R)
        cmd3 = (
            "Rscript " + config.MAT_CORR
            + " " + outdir + "/" + signature_name + "/" + signature_name + ".sorted.txt"
            + " " + outdir + "/" + signature_name
            + "/refinedBySample." + fin_PSI_matrix.split("/")[-1].rstrip(".txt")
            + "." + signature_name + ".sorted.txt"
            + " " + outdir + "/" + signature_name
        )
        print(cmd3)
        os.system(cmd3)

    if GO_analysis_only:
        print("Correlation analysis skipped.")

    print("###GO analysis")
    if GO_analysis:
        getSGlist(outdir, signature_name)
        bg_list = getBGlist(outdir, signature_name)

    if GO_analysis or GO_analysis_only:
        db_versions = args.GOdb_versions
        os.system("mkdir -p " + outdir + "/" + signature_name + "/GO")
        cmd4 = (
            "Rscript " + config.GO_PLOT
            + " " + config.GO_PLOT_LIB
            + " " + bg_list
            + " " + outdir + "/" + signature_name + "/" + signature_name + "_sig_list.txt"
            + " " + outdir + "/" + signature_name + "/GO"
            + " " + signature_name
            + " " + db_versions
        )
        print(cmd4)
        os.system(cmd4)

    if not GO_analysis and not GO_analysis_only:
        print("GO analysis skipped.")


if __name__ == "__main__":
    main()
