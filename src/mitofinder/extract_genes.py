from mitofinder import geneChecker, genbankOutput, runMegahit, circularizationCheck

from Bio import SeqIO, SeqFeature, SeqUtils

from subprocess import Popen
import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys


def main():
    record = SeqIO.read(sys.argv[1], "genbank")
    ID = sys.argv[3]

    with open(sys.argv[2], "a") as importantFeaturesFile:
        for feature in record.features:
            if feature.type.lower() == "cds":
                if "gene" in feature.qualifiers:
                    featureName = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    featureName = feature.qualifiers["product"][0]
                    featureName = "".join(featureName.split())

                importantFeaturesFile.write(">" + str(ID) + "@" + featureName + "\n")
                importantFeaturesFile.write(str(feature.extract(record).seq) + "\n")

            if feature.type.lower() == "rrna":
                if "gene" in feature.qualifiers:
                    featureName = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    featureName = feature.qualifiers["product"][0]
                    featureName = "".join(featureName.split())

                importantFeaturesFile.write(">" + str(ID) + "@" + featureName + "\n")
                importantFeaturesFile.write(str(feature.extract(record).seq) + "\n")


if __name__ == "__main__":
    main()
