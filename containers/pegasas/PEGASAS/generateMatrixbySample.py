import sys
import csv


def loadSampleOrder(fin):
    sample_order = []
    for line in open(fin):
        sample_order = line.strip().split("\t")[1:]
        break
    return sample_order, fin.split("/")[-1].split(".txt")[0]


def readMatrix(fin, sample_order, name, outdir):
    fout = open(
        outdir + "/refinedBySample." + fin.split("/")[-1].split(".txt")[0] + "." + name + ".txt",
        "w",
    )
    header = ""
    for row in open(fin):
        header = row.strip().split("\t")[0:8]
        break
    fout.write("\t".join(header + sample_order) + "\n")
    for row in csv.DictReader(open(fin), dialect="excel-tab"):
        key = "\t".join(
            [row["AC"], row["GeneName"], row["chr"], row["strand"],
             row["exonStart"], row["exonEnd"], row["upstreamEE"], row["downstreamES"]]
        )
        line = key
        for s in sample_order:
            line += "\t" + row[s]
        fout.write(line + "\n")
    fout.close()


def main():
    sample_order, name = loadSampleOrder(sys.argv[2])
    outdir = "/".join(sys.argv[2].split("/")[:-1])
    readMatrix(sys.argv[1], sample_order, name, outdir)


if __name__ == "__main__":
    main()
