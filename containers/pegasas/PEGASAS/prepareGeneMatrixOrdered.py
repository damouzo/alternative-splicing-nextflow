import sys
import os


def loadGeneList(fin):
    """Load ENSEMBL ID / HUGO gene name pairs from CSV."""
    Gene_list = {}
    for line in open(fin):
        ls = line.strip().split(",")
        Gene_list[ls[1]] = ""
    return Gene_list


def loadfromExpMatrix(fin, Gene_list):
    EXP = {}
    n = 0
    header = []
    for line in open(fin):
        if n == 0:
            header = line.strip().split("\t")[1:]
            n += 1
            continue
        ls = line.strip().split("\t")
        l_dict = dict(zip(header, ls[1:]))
        if ls[0].split(".")[0] in Gene_list:
            gene_id = ls[0].split(".")[0]
            if gene_id in EXP:
                print("! Error in", ls[0])
            EXP[gene_id] = {}
            for sample in l_dict.keys():
                EXP[gene_id][sample] = float(l_dict[sample])
    return EXP


def loadOrder(fin):
    order = []
    for line in open(fin):
        order = line.strip().split(",")
        break
    return order


def loadfromGSEA(fin):
    """Load PEGASAS pathway scores; format: sample\tgroup\tscore."""
    ES = {}
    for line in open(fin):
        ls = line.strip().split("\t")
        if ls[1] not in ES:
            ES[ls[1]] = {}
        ES[ls[1]][ls[0]] = float(ls[2])
    return ES


def main():
    fin_gene_score = sys.argv[1]
    Gene_score_dict = loadfromGSEA(fin_gene_score)
    outdir = sys.argv[2].rstrip("/")
    folder_prefix = fin_gene_score.split("/")[-1].split(".scores.txt")[0]
    os.system("mkdir -p " + outdir + " " + outdir + "/" + folder_prefix)
    fout_name = fin_gene_score.split("/")[-1].split(".scores.txt")[0] + ".sorted.txt"
    fout = open(outdir + "/" + folder_prefix + "/" + fout_name, "w")

    header_line = "SampleID"
    value_line  = sys.argv[1].split("/")[-1].split(".")[0]
    sample_order_fin = sys.argv[3]
    order = loadOrder(sample_order_fin)

    if len(order) == 0:
        exit("needs to input the sample order for matrix")

    for dataset in order:
        for sample in sorted(Gene_score_dict[dataset].keys(),
                             key=lambda x: Gene_score_dict[dataset][x]):
            header_line += "\t" + sample
            value_line  += "\t" + str(Gene_score_dict[dataset][sample])

    fout.write(header_line + "\n")
    fout.write(value_line  + "\n")
    fout.close()


if __name__ == "__main__":
    main()
