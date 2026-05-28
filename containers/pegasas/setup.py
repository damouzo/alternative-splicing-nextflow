#!/usr/bin/env python3

from setuptools import setup, find_packages


def main():
    setup(
        name="PEGASAS",
        version="1.1.1",
        description="Pathway Enrichment-Guided Activity Study of Alternative Splicing",
        author="Yang Pan",
        author_email="panyang@ucla.edu",
        packages=["PEGASAS", "PEGASAS.data"],
        scripts=["bin/PEGASAS"],
        include_package_data=True,
        package_data={
            "PEGASAS.data": ["hallmarks50.gmt.txt", "hallmarks50-2.gmt.txt"],
            "PEGASAS": [
                "cor_matrix_direct_perm.R",
                "GO_plot.R",
                "GO_enrichr_plot.R",
            ],
        },
        install_requires=["matplotlib", "numpy", "scipy"],
    )


if __name__ == "__main__":
    main()
