#!/usr/bin/env python3
"""Convert an HMMER3 hmmscan --domtblout file into the tab-separated table that
IsoformSwitchAnalyzeR::analyzePFAM() accepts (its 16-column "new style" layout).

analyzePFAM() does not read domtblout directly: it expects the pfam_scan.pl-style
table (seq id, coordinates, hmm accession/name, scores, clan). Clan is looked up
in Pfam-A.clans.tsv; families without a clan get "No_clan".

Output is empty when the domtblout has no hits, so analyzePFAM() is skipped
downstream instead of choking on a header-only file.
"""

import argparse

HEADER = (
    "# <seq id>\t<alignment start>\t<alignment end>\t<envelope start>"
    "\t<envelope end>\t<hmm acc>\t<hmm name>\t<hmm start>\t<hmm end>"
    "\t<hmm length>\t<bit score>\t<Individual E-value>\t<Conditional E-value>"
    "\t<significant>\t<outcompeted>\t<clan>"
)


def load_clans(path):
    clans = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            acc, clan = fields[0], fields[1].strip()
            clans[acc] = clan if clan else "No_clan"
    return clans


def convert(domtblout, clans):
    rows = []
    with open(domtblout) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # 22 fixed columns followed by a free-text description
            f = line.rstrip("\n").split(None, 22)
            if len(f) < 22:
                continue
            hmm_acc = f[1]
            clan = clans.get(hmm_acc.split(".")[0], "No_clan")
            rows.append("\t".join((
                f[3],   # seq id              <- query name
                f[17],  # alignment start     <- ali from
                f[18],  # alignment end       <- ali to
                f[19],  # envelope start      <- env from
                f[20],  # envelope end        <- env to
                hmm_acc,
                f[0],   # hmm name            <- target name
                f[15],  # hmm start           <- hmm from
                f[16],  # hmm end             <- hmm to
                f[2],   # hmm length          <- tlen
                f[13],  # bit score           <- domain score
                f[12],  # Individual E-value  <- i-Evalue
                f[11],  # Conditional E-value <- c-Evalue
                "1",    # significant
                "0",    # outcompeted
                clan,
            )))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domtblout", required=True)
    ap.add_argument("--clans", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = convert(args.domtblout, load_clans(args.clans))
    with open(args.output, "w") as out:
        if rows:
            out.write(HEADER + "\n")
            out.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
