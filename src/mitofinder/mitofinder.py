from mitofinder._version import __version__
from mitofinder.utils import (
    is_avail,
    is_java_installed,
    find_duplicates,
    check_files,
    get_abs_if_found,
    setup_directory,
)
from mitofinder import (
    circularizationCheck,
    genbankOutput,
    geneChecker,
    runIDBA,
    runMegahit,
    runMetaspades,
)

from Bio import SeqIO, SeqFeature, SeqUtils


from argparse import RawTextHelpFormatter
from datetime import datetime
from shutil import copyfile
from subprocess import Popen
import argparse, os, shlex, shutil, sys
import collections
import glob
import operator
import os.path
import subprocess
import time


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, "".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, "".join(seq))


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith("R|"):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


translation_table = """
Organism genetic code following NCBI table (integer):
1. The Standard Code
2. The Vertebrate Mitochondrial Code
3. The Yeast Mitochondrial Code
4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
5. The Invertebrate Mitochondrial Code
6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
9. The Echinoderm and Flatworm Mitochondrial Code
10. The Euplotid Nuclear Code
11. The Bacterial, Archaeal and Plant Plastid Code
12. The Alternative Yeast Nuclear Code
13. The Ascidian Mitochondrial Code
14. The Alternative Flatworm Mitochondrial Code
16. Chlorophycean Mitochondrial Code
21. Trematode Mitochondrial Code
22. Scenedesmus obliquus Mitochondrial Code
23. Thraustochytrium Mitochondrial Code
24. Pterobranchia Mitochondrial Code
25. Candidate Division SR1 and Gracilibacteria Code
"""
code_dict = {
    1: "1. The Standard Code",
    2: "2. The Vertebrate Mitochondrial Code",
    3: "3. The Yeast Mitochondrial Code",
    4: "4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
    5: "5. The Invertebrate Mitochondrial Code",
    6: "6. The Ciliate, Dasycladacean and Hexamita Nuclear Code",
    9: "9. The Echinoderm and Flatworm Mitochondrial Code",
    10: "10. The Euplotid Nuclear Code",
    11: "11. The Bacterial, Archaeal and Plant Plastid Code",
    12: "12. The Alternative Yeast Nuclear Code",
    13: "13. The Ascidian Mitochondrial Code",
    14: "14. The Alternative Flatworm Mitochondrial Code",
    16: "16. Chlorophycean Mitochondrial Code",
    21: "21. Trematode Mitochondrial Code",
    22: "22. Scenedesmus obliquus Mitochondrial Code",
    23: "23. Thraustochytrium Mitochondrial Code",
    24: "24. Pterobranchia Mitochondrial Code",
    25: "25. Candidate Division SR1 and Gracilibacteria Code",
}


def int_range(value):
    ivalue = int(value)
    if ivalue < 1 or ivalue > 25:
        raise argparse.ArgumentTypeError(
            f"{value} is not in the valid range (1-25). \n {translation_table}"
        )
    return ivalue


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Mitofinder is a pipeline to assemble and annotate mitochondrial DNA from trimmed sequencing reads.",
        prog="mitofinder",
        formatter_class=SmartFormatter,
    )
    parser.add_argument(
        "--megahit",
        help="Use Megahit for assembly. (Default)",
        default=True,
        dest="megahit",
        action="store_true",
    )
    parser.add_argument(
        "--idba",
        help="Use IDBA-UD for assembly. ",
        default=False,
        dest="idba",
        action="store_true",
    )
    parser.add_argument(
        "--metaspades",
        help="Use MetaSPAdes for assembly. ",
        default=False,
        dest="metaspades",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--tRNA-annotation",
        choices=["arwen", "mitfi", "trnascan"],
        help="Choose a tRNA annotation tool to use.Options: ('arwen','mitfi','trnascan') Default = mitfi",
        default="mitfi",
        dest="tRNAannotation",
    )
    parser.add_argument(
        "-j",
        "--seqid",
        help="Job ID to be used throughout the process",
        required=True,
        dest="processName",
    )
    parser.add_argument(
        "-1",
        "--Paired-end1",
        help="File with forward paired-end reads",
        default="",
        dest="PE1",
    )
    parser.add_argument(
        "-2",
        "--Paired-end2",
        help="File with reverse paired-end reads",
        default="",
        dest="PE2",
    )
    parser.add_argument(
        "-s", "--Single-end", help="File with single-end reads", default="", dest="SE"
    )
    parser.add_argument(
        "-a",
        "--assembly",
        help="File with your own assembly",
        default="",
        dest="Assembly",
    )
    parser.add_argument(
        "-m",
        "--max-memory",
        help="max memory to use in Go (MEGAHIT or MetaSPAdes)",
        default="",
        dest="mem",
    )
    parser.add_argument(
        "-l",
        "--length",
        help="Shortest contig length to be used (MEGAHIT). Default = 100",
        type=int,
        default=100,
        dest="shortestContig",
    )
    parser.add_argument(
        "-p",
        "--processors",
        help="Number of threads Mitofinder will use at most.",
        type=int,
        default=4,
        dest="processorsToUse",
    )
    parser.add_argument(
        "-r",
        "--refseq",
        help="Reference mitochondrial genome in GenBank format (.gb).",
        required=True,
        dest="refSeqFile",
    )
    parser.add_argument(
        "-e",
        "--blast-eval",
        help="e-value of blast program used for contig identification and annotation. Default = 0.00001",
        type=float,
        default=0.00001,
        dest="blasteVal",
    )
    parser.add_argument(
        "-n",
        "--nwalk",
        help="Maximum number of codon steps to be tested on each size of the gene to find the start and stop codon during the annotation step. Default = 5 (30 bases)",
        type=int,
        default=5,
        dest="nWalk",
    )
    parser.add_argument(
        "--override",
        help="This option forces MitoFinder to override the previous output directory for the selected assembler.",
        default=False,
        dest="override",
        action="store_true",
    )
    parser.add_argument(
        "--adjust-direction",
        help="This option tells MitoFinder to adjust the direction of selected contig(s) (given the reference).",
        default=False,
        dest="direction",
        action="store_true",
    )
    parser.add_argument(
        "--ignore",
        help="This option tells MitoFinder to ignore the non-standart mitochondrial genes.",
        default=False,
        dest="ignore",
        action="store_true",
    )
    parser.add_argument(
        "--new-genes",
        help="This option tells MitoFinder to try to annotate the non-standard animal mitochondrial genes (e.g. rps3 in fungi). If several references are used, make sure the non-standard genes have the same names in the several references",
        default=False,
        dest="newG",
        action="store_true",
    )
    parser.add_argument(
        "--allow-intron",
        help="This option tells MitoFinder to search for genes with introns. Recommendation : Use it on mitochondrial contigs previously found with MitoFinder without this option.",
        default=False,
        dest="gap",
        action="store_true",
    )
    parser.add_argument(
        "--numt",
        help="This option tells MitoFinder to search for both mitochondrial genes and NUMTs. Recommendation : Use it on nuclear contigs previously found with MitoFinder without this option. ",
        default=False,
        dest="numt",
        action="store_true",
    )
    parser.add_argument(
        "--intron-size",
        help="Size of intron allowed. Default = 5000 bp",
        default=5000,
        type=float,
        dest="intronsize",
    )
    parser.add_argument(
        "--max-contig",
        help="Maximum number of contigs matching to the reference to keep. Default = 0 (unlimited)",
        default=0,
        type=int,
        dest="maxContig",
    )
    parser.add_argument(
        "--cds-merge",
        help="This option tells MitoFinder to not merge the exons in the NT and AA fasta files. ",
        default=True,
        dest="merge",
        action="store_false",
    )
    parser.add_argument(
        "--out-gb",
        help="Do not create annotation output file in GenBank format.",
        default=True,
        dest="genbk",
        action="store_false",
    )
    parser.add_argument(
        "--min-contig-size",
        help="Minimum size of a contig to be considered. Default = 1000",
        default=1000,
        type=float,
        dest="MinContigSize",
    )
    parser.add_argument(
        "--max-contig-size",
        help="Maximum size of a contig to be considered. Default = 25000",
        default=25000,
        type=float,
        dest="MaxContigSize",
    )
    parser.add_argument(
        "--rename-contig",
        help="If set then the contigs matching the reference(s) are renamed. Default: False",
        default=False,
        action="store_true",
        dest="rename",
    )
    parser.add_argument(
        "--blast-identity-nucl",
        help="Nucleotide identity percentage for a hit to be retained. Default = 50",
        default=50,
        type=float,
        dest="blastIdentityNucl",
    )
    parser.add_argument(
        "--blast-identity-prot",
        help="Amino acid identity percentage for a hit to be retained. Default = 40",
        default=40,
        type=float,
        dest="blastIdentityProt",
    )
    parser.add_argument(
        "--blast-size",
        help="Percentage of overlap in blast best hit to be retained. Default = 30",
        default=30,
        type=float,
        dest="aligncutoff",
    )
    parser.add_argument(
        "--circular-size",
        help="Size to consider when checking for circularization. Default = 45",
        default=45,
        type=int,
        dest="circularSize",
    )
    parser.add_argument(
        "--circular-offset",
        help="Offset from start and finish to consider when looking for circularization. Default = 200",
        default=200,
        type=int,
        dest="circularOffSet",
    )
    parser.add_argument(
        "-o",
        "--organism",
        type=int_range,
        choices=range(1, 26),
        help=translation_table,
        required=True,
        dest="organismType",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    parser.add_argument(
        "--example",
        help="Print getting started examples",
        default=False,
        dest="example",
        action="store_true",
    )
    parser.add_argument(
        "--citation",
        help="How to cite MitoFinder",
        default=False,
        dest="citation",
        action="store_true",
    )

    args = parser.parse_args()
    return args


def main():
    args = mainArgs()

    module_dir = os.path.dirname(__file__)
    module_dir = os.path.abspath(module_dir)

    blasteVal = args.blasteVal
    initial_path = os.getcwd()

    # Print example mitofinder cmd and exit
    if args.example == True:
        print(
            """
        # For trimmed paired-end reads:
        mitofinder --megahit -j [seqid] \\\n\t-1 [left_reads.fastq.gz] \\\n\t-2 [right_reads.fastq.gz] \\\n\t-r [genbank_reference.gb] \\\n\t-o [genetic_code] \\\n\t-p [threads] \\\n\t-m [memory]\n\n 
        # For trimmed single-end reads:
        mitofinder --megahit -j [seqid] \\\n\t-s [SE_reads.fastq.gz] \\\n\t-r [genbank_reference.gb] \\\n\t-o [genetic_code] \\\n\t-p [threads] \\\n\t-m [memory]\n\n 
        # For one assembly (one or more contig(s))
        mitofinder -j [seqid] \\\n\t-a [assembly.fasta] \\\n\t-r [genbank_reference.gb] \\\n\t-o [genetic_code] \\\n\t-p [threads] \\\n\t-m [memory]
        """
        )
        exit(0)

    if args.citation == True:
        cite_txt = """ If you use MitoFinder, please cite:
        - Allio, R., Schomaker Bastos, A., Romiguier, J., Prosdocimi, F., Nabholz, B. & Delsuc, F. (2020). MitoFinder: Efficient automated large-scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 20, 892-905. https://doi.org/10.1111/1755-0998.13160
        
        Please also cite the following references depending on the option chosen for the assembly step in MitoFinder:
        - Li, D., Luo, R., Liu, C. M., Leung, C. M., Ting, H. F., Sadakane, K., Yamashita, H. & Lam, T. W. (2016). MEGAHIT v1.0: a fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods, 102(6), 3-11.
        - Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834.
        - Peng, Y., Leung, H. C., Yiu, S. M., & Chin, F. Y. (2012). IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. Bioinformatics, 28(11), 1420-1428.
        
        For tRNAs annotation, depending on the option chosen:
        - MiTFi: Juhling, F., Putz, J., Bernt, M., Donath, A., Middendorf, M., Florentz, C., & Stadler, P. F. (2012). Improved systematic tRNA gene annotation allows new insights into the evolution of mitochondrial tRNA structures and into the mechanisms of mitochondrial genome rearrangements. Nucleic acids research, 40(7), 2833-2845.
        - Laslett, D., & Canback, B. (2008). ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.\n- Chan, P. P., & Lowe, T. M. (2019). tRNAscan-SE: searching for tRNA genes in genomic sequences. In Gene Prediction (pp. 1-14). Humana, New York, NY.
        """
        # Replace long dash with dash
        cite_txt = cite_txt.replace("\u2013", "-")
        print(cite_txt)
        exit(0)

    # Print the selected organism type
    print(f"You have selected translation table: {code_dict[args.organismType]}")

    # Create logfile
    Logfile = args.processName + "_MitoFinder.log"
    Logfile = os.path.join(initial_path, Logfile)
    logfile = open(Logfile, "w")
    logfile.write("Command line: %s" % " ".join(sys.argv) + "\n")
    logfile.write(
        "\nStart time : " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n"
    )
    logfile.write("\nJob name = " + args.processName + "\n")

    start_time = datetime.now()
    print("")
    print("Command line: %s" % " ".join(sys.argv))
    print("")
    print("Now running MitoFinder ...\n")
    print("Start time : " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")

    # Print Job ID
    print("Job ID = " + args.processName)

    # Check java is available if using mitfi
    if args.tRNAannotation == "mitfi":
        if not is_java_installed():
            print(
                "\nERROR: java is not installed/loaded.\nPlease install/load java to run MitoFinder with MiTFi."
            )
            logfile.write(
                "\nERROR: java is not installed/loaded.\nPlease install/load java to run MitoFinder with MiTFi."
            )
            exit(1)

    # Check for either short reads OR a mito assembly.
    if args.PE1 == "" and args.PE2 == "" and args.SE == "" and args.Assembly == "":
        print(
            "\nERROR: Read or assembly files are not specified.\n Please, use -1 -2 -s or -a option to specify input data."
        )
        logfile.write(
            "\nERROR: Read or assembly files are not specified.\n Please, use -1 -2 -s or -a option to specify input data.\n"
        )
        exit(1)

    # Check that short reads are not duplicated
    if (args.PE1 or args.PE2 or args.SE) and not args.Assembly:
        dup_items = find_duplicates([args.PE1, args.PE2, args.SE])

        if dup_items:
            print(
                "Warning: Some short read file names occur more than once!", dup_items
            )
            exit(1)

    # Check that short read files exist if args specified
    check_files([args.PE1, args.PE2, args.SE])

    # Check reference annotation file exists
    if args.refSeqFile:
        print(f"Checking for reference annotation: {args.refSeqFile}")
        args.refSeqFile = get_abs_if_found(args.refSeqFile)

    # Check assembly file exists
    if args.Assembly:
        print(f"Checking for assembly: {args.Assembly}")
        args.Assembly = get_abs_if_found(args.Assembly)

        # If Asm provided but reads or assembler opts also set,
        # warn user.
        if (
            args.megahit
            or args.idba
            or args.metaspades
            or args.PE1
            or args.PE2
            or args.SE
        ) and args.Assembly:
            print("Skipping de novo assembly and using existing assembly.")

        # Set de novo assembly options to false if user provided asm found.
        args.megahit = False
        args.idba = False
        args.metaspades = False

    # Set full path to short reads if filename arg provided
    if args.PE1:
        args.PE1 = get_abs_if_found(args.PE1)
        q1 = "q1=" + args.PE1
    if args.PE2:
        args.PE2 = get_abs_if_found(args.PE2)
        q2 = "q2=" + args.PE2
    if args.SE:
        args.SE = get_abs_if_found(args.SE)
        q3 = "q3=" + args.SE  # Changed from q1 to q3 to prevent SE overwrites PE1

    # Check for both sets of paired-end reads
    if args.PE1:
        if args.PE2:
            T = "PE"
        else:
            print(
                "\nERROR: Only a file with forward paired-end reads was specified.\nPlease specify the file with reverse paired-end reads with -2 option.\nIf you want to use single-end reads, please, use -s option."
            )
            logfile.write(
                "\nERROR: Only a file with forward paired-end reads was specified.\nPlease specify the file with reverse paired-end reads with -2 option.\nIf you want to use single-end reads, please, use -s option.\n"
            )
            exit(1)

    if args.PE2:
        if args.PE1:
            T = "PE"
        else:
            print(
                "\nERROR: Only a file with reverse paired-end reads was specified.\nPlease specify the file with forward paired-end reads with -1 option.\nIf you want to use single reads, please, use -s option."
            )
            logfile.write(
                "\nERROR: Only a file with reverse paired-end reads was specified.\nPlease specify the file with forward paired-end reads with -1 option.\nIf you want to use single reads, please, use -s option.\n"
            )
            exit(1)

    # Populate JOB.input file
    if args.PE1 and args.PE2:
        inputfile = open(args.processName + ".input", "w")
        inputfile.write("type=" + T + "\n" + q1 + "\n" + q2)
        inputfile.close()
        args.inputFile = os.path.join(initial_path, args.processName + ".input")

    if args.SE:
        T = "SE"
        # Metaspades requires PE reads
        if args.metaspades == True:
            print(
                "\nERROR: MetaSPAdes cannot be used for assembly from single-end reads. \nUse Megahit or IDBA-UD.\n"
            )
            logfile.write(
                "\nERROR: MetaSPAdes cannot for assembly from single-end reads. \nUse Megahit or IDBA-UD.\n"
            )
            exit(1)

    # Populate JOB.input file
    if args.SE:
        inputfile = open(args.processName + ".input", "w")
        inputfile.write("type=" + T + "\n" + q3)
        inputfile.close()
        args.inputFile = os.path.join(initial_path, args.processName + ".input")

    args.coveCutOff = 7
    if args.numt == True:
        args.numt = 1
    else:
        args.numt = 0

    if args.gap == True:
        args.gap = 1
    else:
        args.gap = 0

    if args.numt == 1 and args.gap == 1:
        args.nWalk = 0

    # Create output directory
    pathtowork = os.path.join(os.getcwd(), args.processName)
    print("Creating Output directory : " + pathtowork)
    setup_directory(pathtowork)

    logfile.write(
        "\nCreating Output directory : "
        + pathtowork
        + "\nAll results will be written here\n"
    )

    # Check for selected tRNA annotation tool
    # Use version of mifi or arwen that are bundled with MitoFinder if they are requested but not found on PATH
    tRNA = args.tRNAannotation
    if args.tRNAannotation == "trnascan":
        is_avail(["trnascan-SE"])
    elif args.tRNAannotation == "arwen":
        if not is_avail(
            ["arwen"], kill=False
        ):  # Not sure what the executable should be called here
            pathToArwenFolder = os.path.abspath(
                os.path.join(module_dir, "../../bin/arwen/")
            )
            print("Using bundled arwen: {pathToArwenFolder}")
    elif args.tRNAannotation == "mitfi":
        if not is_avail(["mitfi.jar"], kill=False):
            pathToMitfiFolder = os.path.abspath(
                os.path.join(module_dir, "../../bin/mitfi/")
            )
            print("Using bundled mitfi.jar: {pathToMitfiFolder}")

    if not args.Assembly:
        # Check if Megahit is on PATH. Default assembler if no asm provided
        if args.megahit:
            is_avail(["megahit"])
        # Check for IDBA-UD on PATH. Alt assembler.
        if args.idba:
            is_avail(["idba"])
        # Check for MetaSPAdes on PATH. Alt assembler.
        if args.metaspades:
            is_avail(["metaspades.py"])

    # just start the variables for future checking
    firstStep = None  # Megahit
    fourthStep = None  # circularization check
    fifthStep = None  # tRNAscan

    # read the refseq and makeblastdb
    if args.refSeqFile != None:
        if args.refSeqFile[-8:] != ".genbank" and args.refSeqFile[-3:] != ".gb":
            print("Reference mitochondrial genome is not in the expected format")
            print("Provide a file in GenBank format (.gb)")
            print("Aborting")
            logfile.write(
                "Reference mitochondrial genome is not in the expected format\n"
                + "Provide a file in GenBank format (.gb)\n"
                + "Aborting\n"
            )
            exit()

        else:
            gbk_filename = args.refSeqFile
            faa_filename = args.refSeqFile.split("/")[-1].split(".")[0] + ".fasta"
            input_handle = open(gbk_filename, "r").read()
            output_handle = open(pathtowork + "/" + faa_filename, "w")
            translatedGene = open(
                pathtowork + "/translated_genes_for_database.fasta", "w"
            )
            contigdatabase = open(pathtowork + "/contig_id_database.fasta", "w")

            record = input_handle
            importantFeaturesFile = output_handle
            listOfImportantFeatures = {}

            recordCount = record.count("LOCUS   ")
            s = 0
            if record.count("LOCUS   ") > 1:
                c = 0
                for line in range(1, record.count("LOCUS   ") + 1):
                    c += 1
                    refSeq = open(gbk_filename).read()
                    x = refSeq.split("LOCUS   ")[c]
                    tmp = open(args.processName + "_tmp.gb", "w")
                    tmp.write("LOCUS   " + x)
                    tmp.close()
                    record = SeqIO.read(open(args.processName + "_tmp.gb"), "genbank")
                    for feature in record.features:
                        if feature.type.lower() == "cds":
                            if "gene" in feature.qualifiers:
                                featureName = feature.qualifiers["gene"][0]
                            elif "product" in feature.qualifiers:
                                featureName = feature.qualifiers["product"][0]
                            featureName = "".join(featureName.split())
                            featureName = featureName.replace("/", "-")

                            if (
                                not "A" in feature.extract(record).seq.upper()
                                and not "C" in feature.extract(record).seq.upper()
                                and not "G" in feature.extract(record).seq.upper()
                                and not "T" in feature.extract(record).seq.upper()
                            ):
                                print(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                                logfile.write(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                            else:
                                importantFeaturesFile.write(
                                    ">" + record.id + "@" + featureName + "\n"
                                )
                                importantFeaturesFile.write(
                                    str(feature.extract(record).seq) + "\n"
                                )
                                s = 1
                            translatedGene.write(
                                ">" + record.id + "@" + featureName + "\n"
                            )

                            if "translation" in feature.qualifiers:
                                translatedGene.write(
                                    str(feature.qualifiers["translation"][0]) + "\n"
                                )
                            else:
                                translatedGene.write(
                                    str(
                                        feature.extract(record).seq.translate(
                                            table=args.organismType, to_stop=True
                                        )
                                    )
                                    + "\n"
                                )
                                print(
                                    "		WARNING: Reference did not specify a CDS translation for %s. MitoFinder is creating its own from refSeq"
                                    % featureName
                                )
                                logfile.write(
                                    "		WARNING: Reference did not specify a CDS translation for %s. MitoFinder is creating its from refSeq"
                                    % featureName
                                    + "\n"
                                )

                        if feature.type.lower() == "rrna":
                            if "gene" in feature.qualifiers:
                                featureName = feature.qualifiers["gene"][0]
                                featureName = "".join(featureName.split())
                                if (
                                    not "A" in feature.extract(record).seq.upper()
                                    and not "C" in feature.extract(record).seq.upper()
                                    and not "G" in feature.extract(record).seq.upper()
                                    and not "T" in feature.extract(record).seq.upper()
                                ):
                                    print(
                                        "WARNING: no nucleotides sequence have been found for the gene "
                                        + str(featureName)
                                        + " for the reference "
                                        + str(record.id)
                                        + "."
                                    )
                                else:
                                    importantFeaturesFile.write(
                                        ">" + record.id + "@" + featureName + "\n"
                                    )
                                    importantFeaturesFile.write(
                                        str(feature.extract(record).seq) + "\n"
                                    )
                                    listOfImportantFeatures[featureName] = feature
                                    s = 1

                            elif "product" in feature.qualifiers:
                                featureName = feature.qualifiers["product"][0]
                                featureName = "".join(featureName.split())

                                if (
                                    not "A" in feature.extract(record).seq.upper()
                                    and not "C" in feature.extract(record).seq.upper()
                                    and not "G" in feature.extract(record).seq.upper()
                                    and not "T" in feature.extract(record).seq.upper()
                                ):
                                    print(
                                        "WARNING: no nucleotide sequence have been found for the gene "
                                        + str(featureName)
                                        + " for the reference "
                                        + str(record.id)
                                        + "."
                                    )
                                    logfile.write(
                                        "WARNING: no nucleotide sequence have been found for the gene "
                                        + str(featureName)
                                        + " for the reference "
                                        + str(record.id)
                                        + "."
                                    )
                                else:
                                    importantFeaturesFile.write(
                                        ">" + record.id + "@" + featureName + "\n"
                                    )
                                    importantFeaturesFile.write(
                                        str(feature.extract(record).seq) + "\n"
                                    )
                                    listOfImportantFeatures[featureName] = feature
                                    s = 1
                os.remove(args.processName + "_tmp.gb")

            else:
                record = SeqIO.read(open(gbk_filename), "genbank")
                for feature in record.features:
                    if feature.type.lower() == "cds":
                        if "gene" in feature.qualifiers:
                            featureName = feature.qualifiers["gene"][0]
                        elif "product" in feature.qualifiers:
                            featureName = feature.qualifiers["product"][0]
                        featureName = "".join(featureName.split())
                        featureName = featureName.replace("/", "-")
                        if (
                            not "A" in feature.extract(record).seq.upper()
                            and not "C" in feature.extract(record).seq.upper()
                            and not "G" in feature.extract(record).seq.upper()
                            and not "T" in feature.extract(record).seq.upper()
                        ):
                            print(
                                "WARNING: no nucleotide sequence have been found for the gene "
                                + str(featureName)
                                + " for the reference "
                                + str(record.id)
                                + "."
                            )
                            logfile.write(
                                "WARNING: no nucleotide sequence have been found for the gene "
                                + str(featureName)
                                + " for the reference "
                                + str(record.id)
                                + "."
                            )
                        else:
                            importantFeaturesFile.write(
                                ">" + record.id + "@" + featureName + "\n"
                            )
                            importantFeaturesFile.write(
                                str(feature.extract(record).seq) + "\n"
                            )
                            s = 1
                        translatedGene.write(">" + record.id + "@" + featureName + "\n")

                        if "translation" in feature.qualifiers:
                            translatedGene.write(
                                str(feature.qualifiers["translation"][0]) + "\n"
                            )
                        else:
                            translatedGene.write(
                                str(
                                    feature.extract(record).seq.translate(
                                        table=args.organismType, to_stop=True
                                    )
                                )
                                + "\n"
                            )
                            print(
                                "		WARNING: Reference did not specify a CDS translation for %s. MitoFinder is creating its from refSeq"
                                % featureName
                            )
                            logfile.write(
                                "		WARNING: Reference did not specify a CDS translation for %s. MitoFinder is creating its from refSeq"
                                % featureName
                                + "\n"
                            )

                    if feature.type.lower() == "rrna":
                        if "gene" in feature.qualifiers:
                            featureName = feature.qualifiers["gene"][0]
                            featureName = "".join(featureName.split())
                            if (
                                not "A" in feature.extract(record).seq.upper()
                                and not "C" in feature.extract(record).seq.upper()
                                and not "G" in feature.extract(record).seq.upper()
                                and not "T" in feature.extract(record).seq.upper()
                            ):
                                print(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                                logfile.write(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                            else:
                                importantFeaturesFile.write(
                                    ">" + record.id + "@" + featureName + "\n"
                                )
                                importantFeaturesFile.write(
                                    str(feature.extract(record).seq) + "\n"
                                )
                                s = 1
                            listOfImportantFeatures[featureName] = feature

                        elif "product" in feature.qualifiers:
                            featureName = feature.qualifiers["product"][0]
                            featureName = "".join(featureName.split())
                            if (
                                not "A" in feature.extract(record).seq.upper()
                                and not "C" in feature.extract(record).seq.upper()
                                and not "G" in feature.extract(record).seq.upper()
                                and not "T" in feature.extract(record).seq.upper()
                            ):
                                print(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                                logfile.write(
                                    "WARNING: no nucleotide sequence have been found for the gene "
                                    + str(featureName)
                                    + " for the reference "
                                    + str(record.id)
                                    + "."
                                )
                            else:
                                importantFeaturesFile.write(
                                    ">" + record.id + "@" + featureName + "\n"
                                )
                                importantFeaturesFile.write(
                                    str(feature.extract(record).seq) + "\n"
                                )
                                s = 1
                            listOfImportantFeatures[featureName] = feature

            output_handle.close()
            translatedGene.close()
            contigdatabase.close()

            if s == 0:
                print(
                    "\nERROR: MitoFinder didn't found any nucleotide sequence in the reference(s) file(s).\nAborting"
                )
                logfile.write(
                    "ERROR: MitoFinder didn't found any nucleotide sequence in the reference(s) file(s).\nAborting"
                )
                exit()

            # create a gene database
            os.chdir(pathtowork)
            dico_genes = {}
            dico_unknown = {}
            for f in glob.glob("ref*database.fasta*"):
                os.remove(f)
            for f in glob.glob("./*tmp/ref*database.fasta*"):
                os.remove(f)

            for name, seq in read_fasta(open("translated_genes_for_database.fasta")):
                namesp = name
                name = name.split("@")[1]
                if (
                    name.lower() == "coi"
                    or name.lower() == "coxi"
                    or name.lower() == "co1"
                    or name.lower() == "cox1"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubunit1"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubuniti"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubunit1"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubuniti"
                ):
                    name = "COX1"
                elif (
                    name.lower() == "coii"
                    or name.lower() == "coxii"
                    or name.lower() == "co2"
                    or name.lower() == "cox2"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubunit2"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubunitii"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubunit2"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubunitii"
                ):
                    name = "COX2"
                elif (
                    name.lower() == "coiii"
                    or name.lower() == "coxiii"
                    or name.lower() == "co3"
                    or name.lower() == "cox3"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubunit3"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromecoxidasesubunitiii"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubunit3"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "cytochromeoxidasesubunitiii"
                ):
                    name = "COX3"
                elif (
                    name.lower() == "cytb"
                    or name.lower() == "cob"
                    or name.lower().replace(" ", "").replace("-", "") == "cytochromeb"
                ):
                    name = "CYTB"
                elif (
                    name.lower() == "nd1"
                    or name.lower() == "nad1"
                    or name.lower() == "ndh1"
                    or name.lower() == "nadh1"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit1"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubuniti"
                ):
                    name = "ND1"
                elif (
                    name.lower() == "nd2"
                    or name.lower() == "nad2"
                    or name.lower() == "ndh2"
                    or name.lower() == "nadh2"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit2"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitii"
                ):
                    name = "ND2"
                elif (
                    name.lower() == "nd3"
                    or name.lower() == "nad3"
                    or name.lower() == "ndh3"
                    or name.lower() == "nadh3"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit3"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitiii"
                ):
                    name = "ND3"
                elif (
                    name.lower() == "nd4"
                    or name.lower() == "nad4"
                    or name.lower() == "ndh4"
                    or name.lower() == "nadh4"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit4"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitiv"
                ):
                    name = "ND4"
                elif (
                    name.lower() == "nd4l"
                    or name.lower() == "nad4l"
                    or name.lower() == "ndh4l"
                    or name.lower() == "nadh4l"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit4l"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitivl"
                ):
                    name = "ND4L"
                elif (
                    name.lower() == "nd5"
                    or name.lower() == "nad5"
                    or name.lower() == "ndh5"
                    or name.lower() == "nadh5"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit5"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitv"
                ):
                    name = "ND5"
                elif (
                    name.lower() == "nd6"
                    or name.lower() == "nad6"
                    or name.lower() == "ndh6"
                    or name.lower() == "nadh6"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunit6"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "nadhdehydrogenasesubunitvi"
                ):
                    name = "ND6"
                elif (
                    name.lower() == "atp8"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "atpsynthasef0subunit8"
                    or name.lower().replace(" ", "").replace("-", "") == "atpase8"
                ):
                    name = "ATP8"
                elif (
                    name.lower() == "atp6"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "atpsynthasef0subunit6"
                    or name.lower().replace(" ", "").replace("-", "") == "atpase6"
                ):
                    name = "ATP6"
                elif (
                    name.lower().replace(" ", "").replace("-", "") == "rrnl"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "16sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "largesubunitribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "lrrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "16sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "16srrna"
                    or name.lower().replace(" ", "").replace("-", "") == "rnr2"
                    or name.lower().replace(" ", "").replace("-", "") == "mtrnr2"
                    or name.lower().replace(" ", "").replace("-", "") == "rrn16"
                    or name.lower().replace(" ", "").replace("-", "") == "rnl"
                    or name.lower().replace(" ", "").replace("-", "") == "lsu"
                ):
                    name = "rrnL"
                elif (
                    name.lower().replace(" ", "").replace("-", "") == "rrns"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "12sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "smallsubunitribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "srrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "12sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "12srrna"
                    or name.lower().replace(" ", "").replace("-", "") == "rnr1"
                    or name.lower().replace(" ", "").replace("-", "") == "mtrnr1"
                    or name.lower().replace(" ", "").replace("-", "") == "rrn12"
                    or name.lower().replace(" ", "").replace("-", "") == "rns"
                    or name.lower().replace(" ", "").replace("-", "") == "ssu"
                ):
                    name = "rrnS"

                else:
                    if args.ignore == False and args.newG == False:
                        dico_unknown[name] = name
                    elif args.newG == True and not dico_genes.has_key(name):
                        print(
                            'Gene named "'
                            + name
                            + '" in the reference file is not recognized by MitoFinder.'
                        )
                        print("MitoFinder will try to annotate it.")
                        logfile.write(
                            'Gene named "'
                            + name
                            + '" in the reference file is not recognized by MitoFinder.'
                            + "\n"
                            + "MitoFinder will try to annotate it."
                            + "\n\n"
                        )

                    elif args.newG == True and dico_genes.has_key(name):
                        pass
                    elif args.ignore == True and args.newG == False:
                        if not dico_unknown.has_key(name):
                            print(
                                'WARNING: Gene named "'
                                + name
                                + '" in the reference file is not recognized by MitoFinder.'
                            )
                            print("This gene will not be annotated by MitoFinder")
                            logfile.write(
                                'WARNING: Gene named "'
                                + name
                                + '" in the reference file is not recognized by MitoFinder.\n'
                                + "This gene will not be annotated by MitoFinder\n\n"
                            )
                            dico_unknown[name] = name
                        name = "no"

                if name != "rrnL" and name != "rrnS" and name != "no":
                    geneOut = open("ref_" + name + "_database.fasta", "a")
                    geneOut.write(namesp.split("@")[0] + "@" + name + "\n" + seq + "\n")
                    dico_genes[name] = name

            for name, seq in read_fasta(open(faa_filename)):
                namesp = name
                name = name.split("@")[1]
                if (
                    name.lower().replace(" ", "").replace("-", "") == "rrnl"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "16sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "largesubunitribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "lrrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "16sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "16srrna"
                    or name.lower().replace(" ", "").replace("-", "") == "rnr2"
                    or name.lower().replace(" ", "").replace("-", "") == "mtrnr2"
                    or name.lower().replace(" ", "").replace("-", "") == "rrn16"
                    or name.lower().replace(" ", "").replace("-", "") == "rnl"
                    or name.lower().replace(" ", "").replace("-", "") == "lsu"
                ):
                    name = "rrnL"
                elif (
                    name.lower().replace(" ", "").replace("-", "") == "rrns"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "12sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "smallsubunitribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "srrna"
                    or name.lower().replace(" ", "").replace("-", "")
                    == "12sribosomalrna"
                    or name.lower().replace(" ", "").replace("-", "") == "12srrna"
                    or name.lower().replace(" ", "").replace("-", "") == "rnr1"
                    or name.lower().replace(" ", "").replace("-", "") == "mtrnr1"
                    or name.lower().replace(" ", "").replace("-", "") == "rrn12"
                    or name.lower().replace(" ", "").replace("-", "") == "rns"
                    or name.lower().replace(" ", "").replace("-", "") == "ssu"
                ):
                    name = "rrnS"

                else:
                    pass

                if name == "rrnL" or name == "rrnS":
                    geneOut = open("ref_" + name + "_database.fasta", "a")
                    geneOut.write(namesp.split("@")[0] + "@" + name + "\n" + seq + "\n")
                    dico_genes[name] = name

                geneOut = open("contig_id_database.fasta", "a")
                geneOut.write(namesp.split("@")[0] + "@" + name + "\n" + seq + "\n")
                geneOut.close()

            if len(dico_unknown) > 0 and args.ignore == False and args.newG == False:
                if len(dico_unknown) == 1:
                    print(
                        'ERROR: Gene named "'
                        + next(iter(dico_unknown))
                        + '" in the reference file(s) is not recognized by MitoFinder.'
                    )
                    print(
                        "This gene is not a standard mitochondrial gene  (use --ignore or --new-genes options) or please change it to one of the following gene names:"
                    )
                    print(
                        "COX1; COX2; COX3; CYTB; ND1; ND2; ND3; ND4; ND4L; ND5; ND6; ATP8; ATP6; rrnL; rrnS"
                    )
                    print(" ")
                    print(
                        "If you decide to use the new-genes option, please be aware that the names of the new genes in your reference(s) should be homogenized to be considered as unique by MitoFinder (example : Gene1 is not equal to gene1 or gene-1 , use one unique name for all equivalent genes in the different references) "
                    )
                    print(" ")
                    print("Aborting")
                    print(" ")
                    logfile.write(
                        'ERROR: Gene named "'
                        + name
                        + '" in the reference file(s) are not recognized by MitoFinder\n'
                        + "This gene is not a standard mitochondrial gene (use --ignore or --new-genes options) or please change it to one of the following gene names:\n"
                        + "COX1; COX2; COX3; CYTB; ND1; ND2; ND3; ND4; ND4L; ND5; ND6; ATP8; ATP6; rrnL; rrnS\n\nIf you decide to use the new-genes option, please be aware that the names of the new genes in your reference(s) should be homogenized to be considered as unique by MitoFinder (example : Gene1 is not equal to gene1 or gene-1 , use one unique name for all equivalent genes in the different references)\n"
                        + "Aborting\n\n"
                    )
                else:
                    print(
                        "ERROR: The following genes in the reference file(s) are not recognized by MitoFinder."
                    )
                    for k, v in dico_unknown.items():
                        print(" -" + k)
                        logfile.write(" -" + k + "\n")
                    print(
                        "These genes are not standard mitochondrial genes (use --ignore or --new-genes options) or please change them to one of the following gene names:"
                    )
                    print(
                        "COX1; COX2; COX3; CYTB; ND1; ND2; ND3; ND4; ND4L; ND5; ND6; ATP8; ATP6; rrnL; rrnS"
                    )
                    print(" ")
                    print(
                        "If you decide to use the new-genes option, please be aware that the names of the new genes in your reference(s) should be homogenized to be considered as unique by MitoFinder (example : Gene1 is not equal to gene1 or gene-1 , use one unique name for all equivalent genes in the different references) "
                    )
                    print(" ")
                    print("Aborting")
                    print(" ")
                    logfile.write(
                        'ERROR: Gene named "'
                        + name
                        + '" in the reference file(s) are not recognized by MitoFinder\n'
                        + "These genes are not standard mitochondrial genes (use --ignore or --new-genes options) or please change them to one of the following gene names:\n"
                        + "COX1; COX2; COX3; CYTB; ND1; ND2; ND3; ND4; ND4L; ND5; ND6; ATP8; ATP6; rrnL; rrnS\n\nIf you decide to use the new-genes option, please be aware that the names of the new genes in your reference(s) should be homogenized to be considered as unique by MitoFinder (example : Gene1 is not equal to gene1 or gene-1 , use one unique name for all equivalent genes in the different references)\n"
                        + "Aborting\n\n"
                    )
                exit()

            print("")
            logfile.write("\n")

            command = (
                "makeblastdb -in " + str(faa_filename) + " -dbtype nucl"
            )  # need to formatdb refseq first
            args1 = shlex.split(command)
            formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
            formatDB.wait()

            command = "makeblastdb -in contig_id_database.fasta -dbtype nucl"  # need to formatdb refseq first
            args1 = shlex.split(command)
            formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
            formatDB.wait()

    geneList = open(pathtowork + "/genes_list", "w")
    for cle, valeur in dico_genes.items():
        geneList.write(cle + "\n")
    geneList.close()
    for i in (
        "COX1",
        "COX2",
        "COX3",
        "CYTB",
        "ND1",
        "ND2",
        "ND3",
        "ND4",
        "ND4L",
        "ND5",
        "ND6",
        "ATP6",
        "ATP8",
        "rrnL",
        "rrnS",
    ):
        if not i in open(pathtowork + "/genes_list", "r").read():
            print(
                "WARNING: "
                + i
                + " is not in the reference file. MitoFinder will not annotate this gene."
            )
            logfile.write(
                "WARNING: "
                + i
                + " is not in the reference file. MitoFinder will not annotate this gene.\n"
            )

    Assembly = True
    if args.Assembly != "":
        args.megahit = False
        args.idba = False
        args.metaspades = False
        if os.path.exists(pathtowork + "/" + args.processName + "_link.scafSeq"):
            os.remove(pathtowork + "/" + args.processName + "_link.scafSeq")
        os.symlink(args.Assembly, pathtowork + "/" + args.processName + "_link.scafSeq")
        link_file = args.processName + "_link.scafSeq"

    assembler = ""
    if args.megahit == True and args.idba == False and args.metaspades == False:
        assembler = "_megahit"
    if args.idba == True:
        assembler = "_idba"
    if args.metaspades == True:
        assembler = "_metaspades"
    pathOfFinalResults = (
        pathtowork
        + "/"
        + args.processName
        + "_MitoFinder"
        + assembler
        + "_"
        + tRNA
        + "_Final_Results/"
    )
    if os.path.exists(pathOfFinalResults):
        shutil.rmtree(pathOfFinalResults)
    if not os.path.exists(pathOfFinalResults):
        os.makedirs(pathOfFinalResults)

    if Assembly == True:
        logfile.close()
        # let's call megahit
        if args.megahit == True and args.idba == False and args.metaspades == False:
            firstStep = runMegahit.runMegahit(
                processName=args.processName,
                inputFile=args.inputFile,
                shortestContig=args.shortestContig,
                processorsToUse=args.processorsToUse,
                refSeqFile=args.refSeqFile,
                organismType=args.organismType,
                maxMemory=args.mem,
                logfile=Logfile,
                override=args.override,
            )
            out = args.processName + "_megahit"
            logfile = open(Logfile, "a")
            if (
                not os.path.isfile(pathtowork + "/" + out + "/" + out + ".contigs.fa")
                == True
            ):
                print("\n Megahit didn't run")
                print(
                    "Delete or rename the Megahit result folder and restart MitoFinder"
                )
                logfile.write(
                    "\n Megahit didn't run\n"
                    + "Delete or rename the Megahit result folder and restart MitoFinder\n"
                )
                exit()
            if os.path.exists(
                pathtowork + "/" + args.processName + "_link_megahit.scafSeq"
            ):
                os.remove(pathtowork + "/" + args.processName + "_link_megahit.scafSeq")
            os.symlink(
                pathtowork + "/" + out + "/" + out + ".contigs.fa",
                pathtowork + "/" + args.processName + "_link_megahit.scafSeq",
            )
            link_file = args.processName + "_link_megahit.scafSeq"

        # let's call IDBA-UD
        if args.idba == True:
            firstStep = runIDBA.runIDBA(
                processName=args.processName,
                inputFile=args.inputFile,
                shortestContig=args.shortestContig,
                processorsToUse=args.processorsToUse,
                refSeqFile=args.refSeqFile,
                organismType=args.organismType,
                logfile=Logfile,
                override=args.override,
            )
            out = args.processName + "_idba"
            logfile = open(Logfile, "a")
            if not os.path.isfile(pathtowork + "/" + out + "/contig.fa") == True:
                print("\n IDBA-UD didn't run")
                print(
                    "Delete or rename the IDBA-UD result folder and restart MitoFinder"
                )
                logfile.write(
                    "\n IDBA-UD didn't run\n"
                    + "Delete or rename the IDBA-UD result folder and restart MitoFinder\n"
                )
                exit()
            if os.path.exists(
                pathtowork + "/" + args.processName + "_link_idba.scafSeq"
            ):
                os.remove(pathtowork + "/" + args.processName + "_link_idba.scafSeq")
            os.symlink(
                pathtowork + "/" + out + "/contig.fa",
                pathtowork + "/" + args.processName + "_link_idba.scafSeq",
            )
            link_file = args.processName + "_link_idba.scafSeq"

        # let's call MetaSPAdes
        if args.metaspades == True:
            firstStep = runMetaspades.runMetaspades(
                processName=args.processName,
                inputFile=args.inputFile,
                shortestContig=args.shortestContig,
                processorsToUse=args.processorsToUse,
                refSeqFile=args.refSeqFile,
                organismType=args.organismType,
                maxMemory=args.mem,
                logfile=Logfile,
                override=args.override,
            )
            out = args.processName + "_metaspades"
            logfile = open(Logfile, "a")
            if (
                not os.path.isfile(pathtowork + "/" + out + "/" + "scaffolds.fasta")
                == True
            ):
                print("\n MetaSPAdes didn't run")
                print(
                    "Delete or rename the MetaSPAdes result folder and restart MitoFinder"
                )
                logfile.write(
                    "\n MetaSPAdes didn't run\n"
                    + "Delete or rename the MetaSPAdes result folder and restart MitoFinder\n"
                )
                exit()
            if os.path.exists(
                pathtowork + "/" + args.processName + "_link_metaspades.scafSeq"
            ):
                os.remove(
                    pathtowork + "/" + args.processName + "_link_metaspades.scafSeq"
                )
            os.symlink(
                pathtowork + "/" + out + "/" + "scaffolds.fasta",
                pathtowork + "/" + args.processName + "_link_metaspades.scafSeq",
            )
            link_file = args.processName + "_link_metaspades.scafSeq"

        # identification of contigs matching on the refSeq
        blasteVal = args.blasteVal
        if args.metaspades == False and args.megahit == False and args.idba == False:
            logfile = open(Logfile, "a")

        print("Formatting database for mitochondrial contigs identification...")
        logfile.write(
            "Formatting database for mitochondrial contigs identification...\n"
        )
        command = "makeblastdb -in " + link_file + " -dbtype nucl"

        args1 = shlex.split(command)
        formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
        formatDB.wait()

        print("Running mitochondrial contigs identification step...")
        logfile.write("Running mitochondrial contigs identification step...\n")
        with open(args.processName + "_blast_out.txt", "w") as BlastResult:
            command = (
                "blastn -db "
                + link_file
                + " -query contig_id_database.fasta -evalue "
                + str(blasteVal)
                + " -outfmt 6 -perc_identity "
                + str(args.blastIdentityNucl)
            )
            args1 = shlex.split(command)
            blast = Popen(args1, stdout=BlastResult)
            blast.wait()

        os.rename(
            args.processName + "_blast_out.txt",
            pathtowork + "/" + args.processName + "_blast_out.txt",
        )

        mitoblast = open(pathtowork + "/" + args.processName + "_blast_out.txt")

        # fl = sum(1 for line in open(pathtowork+"/"+args.processName+'_blast_out.txt'))
        dico_size_contig = {}

        for r in SeqIO.parse(pathtowork + "/" + link_file, "fasta"):
            dico_size_contig[r.id] = len(r.seq)

        dico_direction = {}
        dico_score = {}
        sup = 0
        for line in mitoblast:
            # testedGene=line.split("\t")[0].split("@")[1]
            contig = line.split("\t")[1].split("\t")[0]
            if (
                float(dico_size_contig.get(contig)) >= args.MinContigSize
                and float(dico_size_contig.get(contig)) <= args.MaxContigSize
            ):
                score = float(line.split("\t")[11])
                ident = line.split("\t")[2]
                start = float(line.split("\t")[8])
                end = float(line.split("\t")[9])
                if start < end:
                    direction = "+"
                else:
                    direction = "-"
                if dico_direction.has_key(contig):
                    dico_direction[contig] = (
                        dico_direction.get(contig) + ";" + direction
                    )
                else:
                    dico_direction[contig] = direction
                if dico_score.has_key(contig):
                    if score > dico_score.get(contig):
                        dico_score[contig] = score
                else:
                    dico_score[contig] = score

            elif (
                float(dico_size_contig.get(contig)) >= args.MinContigSize
                and float(dico_size_contig.get(contig)) >= args.MaxContigSize
            ):
                sup = 1

        sorted_y = sorted(dico_score.items(), key=operator.itemgetter(1), reverse=True)
        sorted_dico_score = collections.OrderedDict(sorted_y)

        if len(collections.OrderedDict(sorted_y)) == 0:
            if sup == 0:
                print("MitoFinder dit not found any mitochondrial contig")
                print("")
                logfile.write("MitoFinder dit not found any mitochondrial contig\n\n")
                time = datetime.now() - start_time
                print("Total wall-clock time used by MitoFinder = " + str(time))
                logfile.write(
                    "\nTotal wall-clock time used by MitoFinder = " + str(time) + "\n"
                )
                exit()
            elif sup == 1:
                print(
                    "MitoFinder dit not found any mitochondrial contig with a size lower than "
                    + str(args.MaxContigSize)
                    + " bp."
                )
                print("")
                logfile.write(
                    "MitoFinder dit not found any mitochondrial contig with a size lower than "
                    + str(args.MaxContigSize)
                    + " bp.\n\n"
                )
                time = datetime.now() - start_time
                print("Total wall-clock time used by MitoFinder = " + str(time))
                logfile.write(
                    "\nTotal wall-clock time used by MitoFinder = " + str(time) + "\n"
                )
                exit()
        elif len(collections.OrderedDict(sorted_y)) == 1:
            print("")
            print("MitoFinder found a single mitochondrial contig")
            print("Checking resulting contig for circularization...")
            logfile.write(
                "\nMitoFinder found a single mitochondrial contig\n"
                + "Checking resulting contig for circularization...\n"
            )
        elif len(collections.OrderedDict(sorted_y)) > 1:
            print("")
            print(
                "MitoFinder found "
                + str(len(collections.OrderedDict(sorted_y)))
                + " contigs matching provided mitochondrial reference(s)"
            )
            if args.maxContig == 0:
                print("Did not check for circularization")
            logfile.write(
                "\nMitoFinder found "
                + str(len(collections.OrderedDict(sorted_y)))
                + " contigs matching provided mitochondrial reference(s)"
                + "\nDid not check for circularization\n"
            )

        if (
            args.maxContig != 0
            and len(collections.OrderedDict(sorted_y)) > args.maxContig
        ):
            if args.maxContig == 1:
                print(
                    "As requested, only the first contig (best match) will be analysed"
                )
                print("Checking resulting contig for circularization...")
                logfile.write(
                    "As requested, only the first contig (best match) will be analysed"
                    + "\nChecking resulting contig for circularization...\n"
                )
            else:
                print(
                    "As requested, only the "
                    + str(args.maxContig)
                    + " first contigs (sorted by hit score) will be analysed"
                )
                print("Did not check for circularization")
                logfile.write(
                    "As requested, only the "
                    + str(args.maxContig)
                    + " first contigs (sorted by hit score) will be analysed"
                    + "\nDid not check for circularization\n"
                )

        fl = 0
        ID_dico = {}
        for k, v in sorted_dico_score.items():
            fl += 1
            if args.maxContig == 0 or fl <= args.maxContig:
                ID_dico[k] = fl
            else:
                fl -= 1
                break

        dico_final_direction = {}
        for key, values in dico_direction.items():
            if ";" in values:
                plus = values.count("+")
                minus = values.count("-")
                if plus > minus:
                    dico_final_direction[key] = "+"
                else:
                    dico_final_direction[key] = "-"
            else:
                dico_final_direction[key] = values

        if fl == 1:
            fout = open(pathtowork + "/" + args.processName + "_contig.fasta", "w")
            for r in SeqIO.parse(pathtowork + "/" + link_file, "fasta"):
                if ID_dico.has_key(r.id):
                    SeqIO.write(r, fout, "fasta")
            fout.close()

            # Resultfile
            pathOfResult = pathtowork + "/" + args.processName + "_contig.fasta"

            # Call module as script
            command = ["python", "-m", "mitofinder.circularizationCheck"]
            # [resultFile, circularSize, circularOffSet]
            circ_args = [pathOfResult, str(args.circularSize), str(args.circularOffSet)]

            # subprocess.run(["python", "-m", "mitofinder.circularizationCheck"] + arguments)
            result = subprocess.run(
                command + circ_args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
            )
            process_output = result.stdout.decode("utf-8") if result.stdout else None
            fourthStep = (
                process_output.rstrip()
                .replace(" ", "")
                .replace("(", "")
                .replace(")", "")
                .split(",")
            )
            # circularizationcheck will return a tuple with (True, start, end)
            print("")
            logfile.write("\n")

            resultFile = args.processName + ".fasta"

            if fourthStep[0] == "True":
                print("Evidences of circularization were found!")
                print(
                    "Sequence is going to be trimmed according to circularization position. \n"
                )
                print("")
                logfile.write(
                    "Evidences of circularization were found!\n"
                    + "Sequence is going to be trimmed according to circularization position. \n\n"
                )
                with open(
                    resultFile, "w"
                ) as outputResult:  # create draft file to be checked and annotated
                    finalResults = SeqIO.read(open(pathOfResult, "rU"), "fasta")
                    finalResults.seq = finalResults.seq[int(fourthStep[2]) :].upper()
                    count = SeqIO.write(
                        finalResults, outputResult, "fasta"
                    )  # trims according to circularization position
            else:
                print(
                    "Evidences of circularization could not be found, but everyother step was successful"
                )
                print("")
                logfile.write(
                    "Evidences of circularization could not be found, but everyother step was successful\n\n"
                )
                with open(
                    resultFile, "w"
                ) as outputResult:  # create draft file to be checked and annotated
                    finalResults = SeqIO.read(open(pathOfResult, "rU"), "fasta")
                    finalResults.seq = finalResults.seq.upper()
                    count = SeqIO.write(
                        finalResults, outputResult, "fasta"
                    )  # no need to trim, since circularization wasn't found

            pathOfFinalResults = (
                pathtowork
                + "/"
                + args.processName
                + "_MitoFinder"
                + assembler
                + "_"
                + tRNA
                + "_Final_Results/"
            )
            if not os.path.exists(pathOfFinalResults):
                os.makedirs(pathOfFinalResults)

            # creating some stat file:
            print("\n")
            print("Creating summary statistics for the mtDNA contig")
            logfile.write("\n\n" + "Creating summary statistics for the mtDNA contig\n")

            finalResults = SeqIO.read(open(resultFile, "rU"), "fasta")
            finalStatsFile = open(pathOfFinalResults + args.processName + ".infos", "w")

            finalStatsFile.write(
                "Initial contig name: " + str(finalResults.id) + "\n\n"
            )

            finalStatsFile.write("Statistics for final sequence:\n")
            finalStatsFile.write("Length: " + str(len(finalResults.seq)) + "\n")
            finalStatsFile.write(
                "GC content: "
                + ("{0:.2f}".format(SeqUtils.GC(finalResults.seq)))
                + "%\n"
            )
            if fourthStep[0] == "True":
                finalStatsFile.write("Circularization: Yes\n")
            else:
                finalStatsFile.write("Circularization: Not found\n")

            if args.direction == True:
                # Call module as script
                command = ["python", "-m", "mitofinder.rename_fasta_seqID"]

                rename_fasta_args = [
                    args.processName,
                    resultFile,
                    os.path.join(
                        pathOfFinalResults, args.processName + "_mtDNA_contig.fasta"
                    ),
                    str(1),
                    dico_final_direction.get(finalResults.id) + str(args.rename),
                ]

            else:
                command = ["python", "-m", "mitofinder.rename_fasta_seqID"]
                rename_fasta_args = (
                    +args.processName
                    + " "
                    + resultFile
                    + " "
                    + pathOfFinalResults
                    + "/"
                    + args.processName
                    + "_mtDNA_contig.fasta"
                    + " "
                    + str(1)
                    + " + "
                    + str(args.rename)
                )
                rename_fasta_args = shlex.split(rename_fasta_args)

            rename = Popen(command + rename_fasta_args, stdout=open(os.devnull, "wb"))
            rename.wait()
            os.remove(resultFile)

            # Annotating with gene_checker
            print("")
            print("Annotating mitochondrial contig")
            print("")
            logfile.write("\nAnnotating\n\n")

            if recordCount > 1:  # if more than 1 ref
                if os.path.isfile(pathtowork + "/ref_for_mtDNA_contig.fasta") == True:
                    os.remove(pathtowork + "/ref_for_mtDNA_contig.fasta")

                for line in open(pathtowork + "/genes_list"):
                    if line.rstrip() != "rrnL" and line.rstrip() != "rrnS":
                        gene = line.rstrip()
                        command = (
                            "makeblastdb -in "
                            + pathtowork
                            + "/ref_"
                            + str(gene + "_database.fasta")
                            + " -dbtype prot"
                        )  # need to formatdb refseq first
                        args1 = shlex.split(command)
                        formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
                        formatDB.wait()
                        with open(
                            pathtowork + "/" + gene + "_blast_out.txt", "w"
                        ) as BlastResultGene:
                            command = (
                                "blastx -db "
                                + "ref_"
                                + gene
                                + "_database.fasta"
                                + " -query "
                                + pathOfFinalResults
                                + "/"
                                + args.processName
                                + "_mtDNA_contig.fasta"
                                + " -evalue "
                                + str(blasteVal)
                                + " -outfmt 6"
                                + " -query_gencode "
                                + str(args.organismType)
                                + " -seg no"
                            )
                            args1 = shlex.split(command)
                            blast = Popen(args1, stdout=BlastResultGene)
                            blast.wait()
                    if line.rstrip() == "rrnL" or line.rstrip() == "rrnS":
                        gene = line.rstrip()
                        command = (
                            "makeblastdb -in "
                            + pathtowork
                            + "/ref_"
                            + str(gene + "_database.fasta")
                            + " -dbtype nucl"
                        )  # need to formatdb refseq first
                        args1 = shlex.split(command)
                        formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
                        formatDB.wait()
                        with open(
                            pathtowork + "/" + gene + "_blast_out.txt", "w"
                        ) as BlastResultGene:
                            command = (
                                "blastn -db "
                                + "ref_"
                                + gene
                                + "_database.fasta"
                                + " -query "
                                + pathOfFinalResults
                                + "/"
                                + args.processName
                                + "_mtDNA_contig.fasta"
                                + " -evalue "
                                + str(blasteVal)
                                + " -outfmt 6 -perc_identity "
                                + str(args.blastIdentityNucl)
                                + " -dust no"
                            )
                            args1 = shlex.split(command)
                            blast = Popen(args1, stdout=BlastResultGene)
                            blast.wait()

                    dico_query = {}
                    bestScore = 0
                    for line in open(pathtowork + "/" + gene + "_blast_out.txt"):
                        query = line.split("\t")[1].split("\t")[0]
                        testedGene = line.split("\t")[1].split("@")[1].split("\t")[0]
                        score = line.split("\t")[11]
                        if dico_query.has_key(testedGene):
                            if float(score) > float(bestScore):
                                dico_query[testedGene] = query
                                bestScore = score
                        else:
                            dico_query[testedGene] = query
                            bestScore = score

                    refFile = open(pathtowork + "/ref_for_mtDNA_contig.fasta", "a")
                    for cle, valeur in dico_query.items():
                        for name, seq in read_fasta(
                            open(pathtowork + "/ref_" + gene + "_database.fasta")
                        ):
                            if name.replace(">", "") == valeur:
                                refFile.write(name + "\n" + seq + "\n")

                    refFile.close()

                if args.gap == 1 or args.numt == 1:
                    command = (
                        "python -m mitofinder.geneChecker_fasta_gaps"
                        + " "
                        + pathtowork
                        + "/ref_for_mtDNA_contig.fasta"
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig.fasta"
                        + " "
                        + args.processName
                        + "_mtDNA_contig.gb"
                        + " "
                        + str(args.organismType)
                        + " "
                        + str(args.aligncutoff)
                        + " "
                        + str(args.coveCutOff)
                        + " "
                        + str(blasteVal)
                        + " "
                        + str(args.blastIdentityProt)
                        + " "
                        + str(args.blastIdentityNucl)
                        + " "
                        + str(args.genbk)
                        + " "
                        + str(args.nWalk)
                        + " "
                        + str(args.intronsize)
                        + " "
                        + str(args.numt)
                        + " "
                        + str(args.gap)
                        + " "
                        + tRNA
                    )
                    args1 = shlex.split(command)
                    fifthStep = Popen(
                        args1,
                        cwd=pathOfFinalResults,
                        stdout=open("geneChecker.log", "a"),
                        stderr=open("geneChecker_error.log", "a"),
                    )
                    fifthStep.wait()
                else:
                    command = (
                        "python -m mitofinder.geneChecker_fasta"
                        + " "
                        + pathtowork
                        + "/ref_for_mtDNA_contig.fasta"
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig.fasta"
                        + " "
                        + args.processName
                        + "_mtDNA_contig.gb"
                        + " "
                        + str(args.organismType)
                        + " "
                        + str(args.aligncutoff)
                        + " "
                        + str(args.coveCutOff)
                        + " "
                        + str(blasteVal)
                        + " "
                        + str(args.blastIdentityProt)
                        + " "
                        + str(args.blastIdentityNucl)
                        + " "
                        + str(args.genbk)
                        + " "
                        + str(args.nWalk)
                        + " "
                        + tRNA
                    )
                    args1 = shlex.split(command)
                    fifthStep = Popen(
                        args1,
                        cwd=pathOfFinalResults,
                        stdout=open("geneChecker.log", "a"),
                        stderr=open("geneChecker_error.log", "a"),
                    )
                    fifthStep.wait()

            else:
                best_ref = open(pathtowork + "/ref_for_mtDNA_contig.fasta", "w")
                for line in open(pathtowork + "/genes_list"):
                    gene = line.rstrip()
                    for name, seq in read_fasta(
                        open(pathtowork + "/ref_" + gene + "_database.fasta")
                    ):
                        best_ref.write(name + "\n" + seq + "\n")
                best_ref.close()

                if args.gap == 1 or args.numt == 1:
                    command = (
                        "python -m mitofinder.geneChecker_fasta_gaps"
                        + " "
                        + pathtowork
                        + "/ref_for_mtDNA_contig.fasta"
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig.fasta"
                        + " "
                        + args.processName
                        + "_mtDNA_contig.gb"
                        + " "
                        + str(args.organismType)
                        + " "
                        + str(args.aligncutoff)
                        + " "
                        + str(args.coveCutOff)
                        + " "
                        + str(blasteVal)
                        + " "
                        + str(args.blastIdentityProt)
                        + " "
                        + str(args.blastIdentityNucl)
                        + " "
                        + str(args.genbk)
                        + " "
                        + str(args.nWalk)
                        + " "
                        + str(args.intronsize)
                        + " "
                        + str(args.numt)
                        + " "
                        + str(args.gap)
                        + " "
                        + tRNA
                    )
                    args1 = shlex.split(command)
                    fifthStep = Popen(
                        args1,
                        cwd=pathOfFinalResults,
                        stdout=open("geneChecker.log", "a"),
                        stderr=open("geneChecker_error.log", "a"),
                    )
                    fifthStep.wait()
                else:
                    command = (
                        "python -m mitofinder.geneChecker_fasta"
                        + " "
                        + pathtowork
                        + "/ref_for_mtDNA_contig.fasta"
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig.fasta"
                        + " "
                        + args.processName
                        + "_mtDNA_contig.gb"
                        + " "
                        + str(args.organismType)
                        + " "
                        + str(args.aligncutoff)
                        + " "
                        + str(args.coveCutOff)
                        + " "
                        + str(blasteVal)
                        + " "
                        + str(args.blastIdentityProt)
                        + " "
                        + str(args.blastIdentityNucl)
                        + " "
                        + str(args.genbk)
                        + " "
                        + str(args.nWalk)
                        + " "
                        + tRNA
                    )
                    args1 = shlex.split(command)
                    fifthStep = Popen(
                        args1,
                        cwd=pathOfFinalResults,
                        stdout=open("geneChecker.log", "w"),
                        stderr=open("geneChecker_error.log", "w"),
                    )
                    fifthStep.wait()
            # Check if correct output was created
            if tRNA == "arwen":
                test_arwen = (
                    pathOfFinalResults + "/" + args.processName + "_mtDNA_contig.arwen"
                )
                if os.path.isfile(test_arwen) == True:
                    print("tRNA annotation with Arwen run well.\n")
                    logfile.write("tRNA annotation with Arwen run well.\n\n")
                else:
                    print(
                        "ERROR: tRNA annotation failed.\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n"
                    )
                    logfile.write(
                        "ERROR: tRNA annotation failed.\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n\n"
                    )
                    exit()
            elif tRNA == "trnascan":
                test_trnascan = (
                    pathOfFinalResults
                    + "/"
                    + args.processName
                    + "_mtDNA_contig.trnascan"
                )
                if os.path.isfile(test_trnascan) == True:
                    print("tRNA annotation with tRNAscan-SE run well.\n")
                    logfile.write("tRNA annotation with tRNAscan-SE run well.\n\n")
                else:
                    print(
                        "ERROR: tRNA annotation failed.\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n"
                    )
                    logfile.write(
                        "ERROR: tRNA annotation failed.\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n\n"
                    )
                    exit()
            elif tRNA == "mitfi":

                def check_if_string_in_file(f, string):
                    with open(f, "r") as read_obj:
                        for line in read_obj:
                            if string in line:
                                return True
                    return False

                if check_if_string_in_file(
                    pathtowork + "/geneChecker.log", "MiTFi failed."
                ) or (
                    os.stat(pathOfFinalResults + "MiTFi.log").st_size != 0
                    and not check_if_string_in_file(
                        pathOfFinalResults + "MiTFi.log", "hits"
                    )
                ):
                    print(
                        "ERROR: tRNA annotation failed.\nTo see what happened, please check:\n"
                        + pathtowork
                        + "/geneChecker_error.log \n"
                        + pathtowork
                        + "/geneChecker.log \nor "
                        + pathOfFinalResults
                        + "MiTFi.log\nAborting\n"
                    )
                    logfile.write(
                        "ERROR: tRNA annotation failed.\nTo see what happened, please check:\n"
                        + pathtowork
                        + "/geneChecker_error.log \n"
                        + pathtowork
                        + "/geneChecker.log \nor "
                        + pathOfFinalResults
                        + "MiTFi.log\nAborting\n\n"
                    )
                    exit()
                else:
                    print("tRNA annotation with MitFi run well.\n")
                    logfile.write("tRNA annotation with MitFi run well.\n\n")

            test_gene_checker = (
                pathOfFinalResults + args.processName + "_mtDNA_contig.gb"
            )
            if os.path.isfile(test_gene_checker) == True:
                print("Annotation completed\n")
                logfile.write("Annotation completed\n\n")
            else:
                print(
                    "ERROR: Gene annotation failed\nPlease check  "
                    + pathtowork
                    + "/geneChecker_error.log to see what happened\nAborting\n"
                )
                logfile.write(
                    "ERROR: Gene annotation failed\nPlease check  "
                    + pathtowork
                    + "/geneChecker_error.log to see what happened\nAborting\n\n"
                )
                exit()

        elif fl > 1:
            if os.path.isfile(pathtowork + "/" + "geneChecker.log") == True:
                os.remove(pathtowork + "/" + "geneChecker.log")
            if os.path.isfile(pathtowork + "/" + "geneChecker_error.log") == True:
                os.remove(pathtowork + "/" + "geneChecker_error.log")

            # Extract every contigs one by one
            contg_list = open(pathtowork + "/" + "contig_list.txt", "w")
            for r in SeqIO.parse(pathtowork + "/" + link_file, "fasta"):
                if ID_dico.has_key(r.id):
                    fout = open(
                        pathtowork
                        + "/"
                        + args.processName
                        + "_contig_"
                        + str(ID_dico.get(r.id))
                        + ".fasta",
                        "w",
                    )
                    contg_list.write(
                        args.processName
                        + "_contig_"
                        + str(ID_dico.get(r.id))
                        + ".fasta"
                        + "\n"
                    )
                    SeqIO.write(r, fout, "fasta")
            fout.close()
            contg_list.close()

            c = 1
            for line in open(pathtowork + "/" + "contig_list.txt", "r"):
                pathOfResult = (
                    pathtowork + "/" + args.processName + "_contig_" + str(c) + ".fasta"
                )

                resultFile = args.processName + "_mtDNA_contig_" + str(c) + ".fasta"

                with open(
                    resultFile, "w"
                ) as outputResult:  # create draft file to be checked and annotated
                    finalResults = SeqIO.read(open(pathOfResult, "rU"), "fasta")
                    finalResults.seq = finalResults.seq.upper()
                    count = SeqIO.write(
                        finalResults, outputResult, "fasta"
                    )  # no need to trim, since circularization wasn't found

                pathOfFinalResults = (
                    pathtowork
                    + "/"
                    + args.processName
                    + "_MitoFinder"
                    + assembler
                    + "_"
                    + tRNA
                    + "_Final_Results/"
                )
                if not os.path.exists(pathOfFinalResults):
                    os.makedirs(pathOfFinalResults)

                # creating some stat file:
                print("\n")
                print("Creating summary statistics for mtDNA contig " + str(c))
                print("")
                logfile.write(
                    "\n\n"
                    + "Creating summary statistics for mtDNA contig "
                    + str(c)
                    + "\n\n"
                )

                finalResults = SeqIO.read(open(resultFile, "rU"), "fasta")
                finalStatsFile = open(
                    pathOfFinalResults
                    + args.processName
                    + "_mtDNA_contig_"
                    + str(c)
                    + ".infos",
                    "w",
                )

                finalStatsFile.write("Statistics for contig " + str(c) + ":\n\n")
                finalStatsFile.write(
                    "Initial contig name: " + str(finalResults.id) + "\n"
                )
                finalStatsFile.write("Length: " + str(len(finalResults.seq)) + "\n")
                finalStatsFile.write(
                    "GC content: "
                    + ("{0:.2f}".format(SeqUtils.GC(finalResults.seq)))
                    + "%\n"
                )
                finalStatsFile.write("Circularization: NA\n")

                if args.direction == True:
                    command = (
                        "python -m mitofinder.rename_fasta_seqID"
                        + " "
                        + args.processName
                        + " "
                        + resultFile
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig_"
                        + str(c)
                        + ".fasta"
                        + " "
                        + str(c)
                        + " "
                        + dico_final_direction.get(finalResults.id)
                        + " "
                        + str(args.rename)
                    )
                else:
                    command = (
                        "python -m mitofinder.rename_fasta_seqID"
                        + " "
                        + args.processName
                        + " "
                        + resultFile
                        + " "
                        + pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig_"
                        + str(c)
                        + ".fasta"
                        + " "
                        + str(c)
                        + " + "
                        + str(args.rename)
                    )
                args1 = shlex.split(command)
                rename = Popen(args1, stdout=open(os.devnull, "wb"))
                rename.wait()

                os.remove(resultFile)
                # shutil.copyfile(resultFile, pathOfFinalResults+"/"+args.processName+"_mtDNA_contig_"+str(c)+".fasta")

                # creating best ref file for annotation

                command = (
                    "makeblastdb -in "
                    + pathOfFinalResults
                    + "/"
                    + args.processName
                    + "_mtDNA_contig_"
                    + str(c)
                    + ".fasta"
                    + " -dbtype nucl"
                )  # need to formatdb refseq first
                args1 = shlex.split(command)
                formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
                formatDB.wait()

                if recordCount > 1:  # if more than 1 ref
                    print("Looking for best reference genes for mtDNA contig " + str(c))
                    print("")
                    logfile.write(
                        "Looking for best reference genes for mtDNA contig "
                        + str(c)
                        + "\n\n"
                    )

                    if (
                        os.path.isfile(
                            pathtowork + "/ref_for_contig_" + str(c) + ".fasta"
                        )
                        == True
                    ):
                        os.remove(pathtowork + "/ref_for_contig_" + str(c) + ".fasta")

                    for line in open(pathtowork + "/genes_list"):
                        if line.rstrip() != "rrnL" and line.rstrip() != "rrnS":
                            gene = line.rstrip()
                            command = (
                                "makeblastdb -in "
                                + pathtowork
                                + "/ref_"
                                + str(gene + "_database.fasta")
                                + " -dbtype prot"
                            )  # need to formatdb refseq first
                            args1 = shlex.split(command)
                            formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
                            formatDB.wait()
                            with open(
                                pathtowork + "/" + gene + "_blast_out.txt", "w"
                            ) as BlastResultGene:
                                command = (
                                    "blastx -db "
                                    + "ref_"
                                    + gene
                                    + "_database.fasta"
                                    + " -query "
                                    + pathOfFinalResults
                                    + "/"
                                    + args.processName
                                    + "_mtDNA_contig_"
                                    + str(c)
                                    + ".fasta"
                                    + " -evalue "
                                    + str(blasteVal)
                                    + " -outfmt 6"
                                    + " -query_gencode "
                                    + str(args.organismType)
                                    + " -seg no"
                                )
                                args1 = shlex.split(command)
                                blast = Popen(args1, stdout=BlastResultGene)
                                blast.wait()
                        if line.rstrip() == "rrnL" or line.rstrip() == "rrnS":
                            gene = line.rstrip()
                            command = (
                                "makeblastdb -in "
                                + pathtowork
                                + "/ref_"
                                + str(gene + "_database.fasta")
                                + " -dbtype nucl"
                            )  # need to formatdb refseq first
                            args1 = shlex.split(command)
                            formatDB = Popen(args1, stdout=open(os.devnull, "wb"))
                            formatDB.wait()
                            with open(
                                pathtowork + "/" + gene + "_blast_out.txt", "w"
                            ) as BlastResultGene:
                                command = (
                                    "blastn -db "
                                    + "ref_"
                                    + gene
                                    + "_database.fasta"
                                    + " -query "
                                    + pathOfFinalResults
                                    + "/"
                                    + args.processName
                                    + "_mtDNA_contig_"
                                    + str(c)
                                    + ".fasta"
                                    + " -evalue "
                                    + str(blasteVal)
                                    + " -outfmt 6 -perc_identity "
                                    + str(args.blastIdentityNucl)
                                    + " -dust no"
                                )
                                args1 = shlex.split(command)
                                blast = Popen(args1, stdout=BlastResultGene)
                                blast.wait()

                        dico_query = {}
                        bestScore = 0
                        for line in open(pathtowork + "/" + gene + "_blast_out.txt"):
                            query = line.split("\t")[1].split("\t")[0]
                            testedGene = (
                                line.split("\t")[1].split("@")[1].split("\t")[0]
                            )
                            score = line.split("\t")[-1]
                            if dico_query.has_key(testedGene):
                                if float(score) > float(bestScore):
                                    dico_query[testedGene] = query
                                    bestScore = score
                            else:
                                dico_query[testedGene] = query
                                bestScore = score

                        refFile = open(
                            pathtowork + "/ref_for_contig_" + str(c) + ".fasta", "a"
                        )
                        for cle, valeur in dico_query.items():
                            for name, seq in read_fasta(
                                open(pathtowork + "/ref_" + gene + "_database.fasta")
                            ):
                                if name.replace(">", "") == valeur:
                                    refFile.write(name + "\n" + seq + "\n")

                        refFile.close()

                print("Annotating mtDNA contig " + str(c))
                print("")
                logfile.write("Annotating mtDNA contig " + str(c) + "\n\n")

                # Annotating with gene_checker

                if recordCount > 1:
                    if args.gap == 1 or args.numt == 1:
                        command = (
                            "python -m mitofinder.geneChecker_fasta_gaps"
                            + " "
                            + pathtowork
                            + "/ref_for_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + pathOfFinalResults
                            + "/"
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".gb"
                            + " "
                            + str(args.organismType)
                            + " "
                            + str(args.aligncutoff)
                            + " "
                            + str(args.coveCutOff)
                            + " "
                            + str(blasteVal)
                            + " "
                            + str(args.blastIdentityProt)
                            + " "
                            + str(args.blastIdentityNucl)
                            + " "
                            + str(args.genbk)
                            + " "
                            + str(args.nWalk)
                            + " "
                            + str(args.intronsize)
                            + " "
                            + str(args.numt)
                            + " "
                            + str(args.gap)
                            + " "
                            + tRNA
                        )
                        args1 = shlex.split(command)
                        fifthStep = Popen(
                            args1,
                            cwd=pathOfFinalResults,
                            stdout=open("geneChecker.log", "a"),
                            stderr=open("geneChecker_error.log", "a"),
                        )
                        fifthStep.wait()
                    else:
                        command = (
                            "python -m mitofinder.geneChecker_fasta"
                            + " "
                            + pathtowork
                            + "/ref_for_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + pathOfFinalResults
                            + "/"
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".gb"
                            + " "
                            + str(args.organismType)
                            + " "
                            + str(args.aligncutoff)
                            + " "
                            + str(args.coveCutOff)
                            + " "
                            + str(blasteVal)
                            + " "
                            + str(args.blastIdentityProt)
                            + " "
                            + str(args.blastIdentityNucl)
                            + " "
                            + str(args.genbk)
                            + " "
                            + str(args.nWalk)
                            + " "
                            + tRNA
                        )
                        args1 = shlex.split(command)
                        fifthStep = Popen(
                            args1,
                            cwd=pathOfFinalResults,
                            stdout=open("geneChecker.log", "a"),
                            stderr=open("geneChecker_error.log", "a"),
                        )
                        fifthStep.wait()

                else:
                    # creation du fichier .fasta pour geneChecker_fasta.py
                    best_ref = open(pathtowork + "/ref_for_contigs.fasta", "w")
                    for line in open(pathtowork + "/genes_list"):
                        gene = line.rstrip()
                        for name, seq in read_fasta(
                            open(pathtowork + "/ref_" + gene + "_database.fasta")
                        ):
                            best_ref.write(name + "\n" + seq + "\n")
                    best_ref.close()

                    if args.gap == 1 or args.numt == 1:
                        command = (
                            "python -m mitofinder.geneChecker_fasta_gaps"
                            + " "
                            + pathtowork
                            + "/ref_for_contigs.fasta"
                            + " "
                            + pathOfFinalResults
                            + "/"
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".gb"
                            + " "
                            + str(args.organismType)
                            + " "
                            + str(args.aligncutoff)
                            + " "
                            + str(args.coveCutOff)
                            + " "
                            + str(blasteVal)
                            + " "
                            + str(args.blastIdentityProt)
                            + " "
                            + str(args.blastIdentityNucl)
                            + " "
                            + str(args.genbk)
                            + " "
                            + str(args.nWalk)
                            + " "
                            + str(args.intronsize)
                            + " "
                            + str(args.numt)
                            + " "
                            + str(args.gap)
                            + " "
                            + tRNA
                        )
                        args1 = shlex.split(command)
                        fifthStep = Popen(
                            args1,
                            cwd=pathOfFinalResults,
                            stdout=open("geneChecker.log", "a"),
                            stderr=open("geneChecker_error.log", "a"),
                        )
                        fifthStep.wait()
                    else:
                        command = (
                            "python -m mitofinder.geneChecker_fasta"
                            + " "
                            + pathtowork
                            + "/ref_for_contigs.fasta"
                            + " "
                            + pathOfFinalResults
                            + "/"
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".fasta"
                            + " "
                            + args.processName
                            + "_mtDNA_contig_"
                            + str(c)
                            + ".gb"
                            + " "
                            + str(args.organismType)
                            + " "
                            + str(args.aligncutoff)
                            + " "
                            + str(args.coveCutOff)
                            + " "
                            + str(blasteVal)
                            + " "
                            + str(args.blastIdentityProt)
                            + " "
                            + str(args.blastIdentityNucl)
                            + " "
                            + str(args.genbk)
                            + " "
                            + str(args.nWalk)
                            + " "
                            + tRNA
                        )
                        args1 = shlex.split(command)
                        fifthStep = Popen(
                            args1,
                            cwd=pathOfFinalResults,
                            stdout=open("geneChecker.log", "w"),
                            stderr=open("geneChecker_error.log", "w"),
                        )
                        fifthStep.wait()

                if tRNA == "arwen":
                    test_arwen = (
                        pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig_"
                        + str(c)
                        + ".arwen"
                    )
                    if os.path.isfile(test_arwen) == True:
                        print("tRNA annotation with Arwen run well.\n")
                        logfile.write("tRNA annotation with Arwen run well.\n\n")
                    else:
                        print(
                            "ERROR: tRNA annotation failed.\nPlease check  "
                            + pathtowork
                            + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n"
                        )
                        logfile.write(
                            "ERROR: tRNA annotation failed.\nPlease check  "
                            + pathtowork
                            + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n\n"
                        )
                        exit()
                elif tRNA == "trnascan":
                    test_trnascan = (
                        pathOfFinalResults
                        + "/"
                        + args.processName
                        + "_mtDNA_contig_"
                        + str(c)
                        + ".trnascan"
                    )
                    if os.path.isfile(test_trnascan) == True:
                        print("tRNA annotation with tRNAscan-SE run well.\n")
                        logfile.write("tRNA annotation with tRNAscan-SE run well.\n\n")
                    else:
                        print(
                            "ERROR: tRNA annotation failed.\nPlease check  "
                            + pathtowork
                            + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n"
                        )
                        logfile.write(
                            "ERROR: tRNA annotation failed.\nPlease check  "
                            + pathtowork
                            + "/geneChecker_error.log or geneChecker.log to see what happened\nAborting\n\n"
                        )
                        exit()
                elif tRNA == "mitfi":

                    def check_if_string_in_file(f, string):
                        with open(f, "r") as read_obj:
                            for line in read_obj:
                                if string in line:
                                    return True
                        return False

                    if check_if_string_in_file(
                        pathtowork + "/geneChecker.log", "MiTFi failed."
                    ) or (
                        os.stat(pathOfFinalResults + "MiTFi.log").st_size != 0
                        and not check_if_string_in_file(
                            pathOfFinalResults + "MiTFi.log", "hits"
                        )
                    ):
                        print(
                            "ERROR: tRNA annotation failed.\nTo see what happened, please check:\n"
                            + pathtowork
                            + "/geneChecker_error.log \n"
                            + pathtowork
                            + "/geneChecker.log \nor "
                            + pathOfFinalResults
                            + "MiTFi.log\nAborting\n"
                        )
                        logfile.write(
                            "ERROR: tRNA annotation failed.\nTo see what happened, please check:\n"
                            + pathtowork
                            + "/geneChecker_error.log \n"
                            + pathtowork
                            + "/geneChecker.log \nor "
                            + pathOfFinalResults
                            + "MiTFi.log\nto see what happened\nAborting\n\n"
                        )
                        exit()
                    else:
                        print("tRNA annotation with MitFi run well.\n")
                        logfile.write("tRNA annotation with MitFi run well.\n\n")
                test_gene_checker = (
                    pathOfFinalResults
                    + args.processName
                    + "_mtDNA_contig_"
                    + str(c)
                    + ".gb"
                )
                if os.path.isfile(test_gene_checker) == True:
                    print("Annotation completed\n")
                    logfile.write("Annotation completed\n" + "\n")
                else:
                    print(
                        "ERROR: Gene annotation failed for mtDNA contig "
                        + str(c)
                        + ".\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log to see what happened\nAborting\n"
                    )
                    logfile.write(
                        "ERROR: Gene annotation failed for mtDNA contig "
                        + str(c)
                        + ".\nPlease check  "
                        + pathtowork
                        + "/geneChecker_error.log to see what happened\nAborting\n"
                        + "\n"
                    )
                    exit()
                c += 1

    # Creating GFF and fasta file

    print("\nCreating GFF and fasta files.\n")
    print("Note: ")
    logfile.write("\nCreating GFF and fasta files.\n\n" + "Note: " + "\n")
    for f in sorted(glob.glob(pathOfFinalResults + "/*.gb"), key=os.path.getmtime):
        gnb = 0
        out_fasta_nt = f.split(".gb")[0] + "_genes_NT.fasta"
        out_fasta_aa = f.split(".gb")[0] + "_genes_AA.fasta"
        out_fasta_nt = open(out_fasta_nt, "w")
        out_fasta_aa = open(out_fasta_aa, "w")
        dgen = {}
        with open(f) as infile:
            record = SeqIO.read(infile, "genbank")
            for feature in record.features:
                seq_aa = ""
                if feature.type.lower() == "cds" or feature.type == "rRNA":
                    if "gene" in feature.qualifiers:
                        featureName = feature.qualifiers["gene"][0]
                        if "translation" in feature.qualifiers:
                            seq_aa = feature.qualifiers["translation"][0]
                    elif "product" in feature.qualifiers:
                        featureName = feature.qualifiers["product"][0]
                        if "translation" in feature.qualifiers:
                            seq_aa = feature.qualifiers["translation"][0]
                    featureName = "".join(featureName.split()).split("_")[0]
                    if dgen.has_key(featureName):
                        dgen[featureName] += 1
                    else:
                        dgen[featureName] = 1
        dico_cds_n = {}
        dico_cds_aa = {}
        dico_cds_nt = {}
        with open(f) as infile:
            record = SeqIO.read(infile, "genbank")
            for feature in record.features:
                seq_aa = ""
                if feature.type.lower() == "cds" or feature.type == "rRNA":
                    if "gene" in feature.qualifiers:
                        featureName = feature.qualifiers["gene"][0]
                        if "translation" in feature.qualifiers:
                            seq_aa = feature.qualifiers["translation"][0]
                    elif "product" in feature.qualifiers:
                        featureName = feature.qualifiers["product"][0]
                        if "translation" in feature.qualifiers:
                            seq_aa = feature.qualifiers["translation"][0]
                    featureName = "".join(featureName.split()).split("_")[0]
                    if dgen.get(featureName) == 1:
                        gnb = gnb + 1
                        out_fasta_nt.write(
                            ">" + args.processName + "@" + featureName + "\n"
                        )
                        out_fasta_nt.write(str(feature.extract(record).seq) + "\n")
                        if str(seq_aa) != "":
                            out_fasta_aa.write(
                                ">" + args.processName + "@" + featureName + "\n"
                            )
                            out_fasta_aa.write(str(seq_aa) + "\n")

                    elif not dico_cds_n.has_key(featureName) and args.merge == True:
                        gnb = gnb + 1
                        dico_cds_n[featureName] = 1
                        dico_cds_nt[featureName] = str(feature.extract(record).seq)
                        if str(seq_aa) != "":
                            dico_cds_aa[featureName] = str(seq_aa)
                    elif (
                        dico_cds_n.has_key(featureName)
                        and dico_cds_n.get(featureName) < dgen.get(featureName)
                        and args.merge == True
                    ):
                        dico_cds_n[featureName] += 1
                        if "(+)" in str(feature.location):
                            dico_cds_nt[featureName] = str(
                                dico_cds_nt.get(featureName)
                            ) + str(feature.extract(record).seq)
                            if str(seq_aa) != "":
                                dico_cds_aa[featureName] = str(
                                    dico_cds_aa.get(featureName)
                                ) + str(seq_aa)
                        else:
                            dico_cds_nt[featureName] = str(
                                feature.extract(record).seq
                            ) + str(dico_cds_nt.get(featureName))
                            if str(seq_aa) != "":
                                dico_cds_aa[featureName] = str(seq_aa) + str(
                                    dico_cds_aa.get(featureName)
                                )
                    else:
                        if not dico_cds_n.has_key(featureName):
                            dico_cds_n[featureName] = 1
                            gnb = gnb + 1
                        out_fasta_nt.write(
                            ">" + args.processName + "@" + featureName + "\n"
                        )
                        out_fasta_nt.write(str(feature.extract(record).seq) + "\n")
                        if str(seq_aa) != "":
                            out_fasta_aa.write(
                                ">" + args.processName + "@" + featureName + "\n"
                            )
                            out_fasta_aa.write(str(seq_aa) + "\n")

                    if (
                        dico_cds_n.has_key(featureName)
                        and dico_cds_n.get(featureName) == dgen.get(featureName)
                        and args.merge == True
                    ):
                        out_fasta_nt.write(
                            ">" + args.processName + "@" + featureName + "\n"
                        )
                        out_fasta_nt.write(str(dico_cds_nt.get(featureName)) + "\n")
                        if str(dico_cds_aa.get(featureName)) != "":
                            out_fasta_aa.write(
                                ">" + args.processName + "@" + featureName + "\n"
                            )
                            out_fasta_aa.write(str(dico_cds_aa.get(featureName)) + "\n")

        out_fasta_nt.close()
        out_fasta_aa.close()
        if gnb == 0:
            print(
                str(gnb)
                + " gene was found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
            )
            logfile.write(
                str(gnb)
                + " gene was found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
                + "\n"
            )
            os.remove(f.split(".gb")[0] + "_genes_NT.fasta")
            os.remove(f.split(".gb")[0] + "_genes_AA.fasta")
        elif gnb == 1:
            print(
                str(gnb)
                + " gene was found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
            )
            logfile.write(
                str(gnb)
                + " gene was found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
                + "\n"
            )
        else:
            print(
                str(gnb)
                + " genes were found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
            )
            logfile.write(
                str(gnb)
                + " genes were found in "
                + f.split(".gb")[0].split("Final_Results/" + args.processName + "_")[1]
                + "\n"
            )

    print("")
    logfile.write("\n")

    # sort gff
    pathlist = glob.glob(pathOfFinalResults + "/*.gb")
    if len(pathlist) == 1:
        for f in sorted(
            glob.glob(pathOfFinalResults + "/*_raw.gff"), key=os.path.getmtime
        ):
            command = (
                "python -m mitofinder.sort_gff"
                + " "
                + f
                + " "
                + args.processName
                + ".1 "
                + str(args.organismType)
                + " "
                + str(args.rename)
            )
            args1 = shlex.split(command)
            sort_gff = Popen(
                args1,
                cwd=pathOfFinalResults,
                stdout=open("geneChecker.log", "a"),
                stderr=open("geneChecker_error.log", "a"),
            )
            sort_gff.wait()
    else:
        for f in sorted(
            glob.glob(pathOfFinalResults + "/*_raw.gff"), key=os.path.getmtime
        ):
            command = (
                "python -m mitofinder.sort_gff"
                + " "
                + f
                + " "
                + args.processName
                + "."
                + f.split("_raw")[0][-1]
                + " "
                + str(args.organismType)
                + " "
                + str(args.rename)
            )
            args1 = shlex.split(command)
            sort_gff = Popen(
                args1,
                cwd=pathOfFinalResults,
                stdout=open("geneChecker.log", "a"),
                stderr=open("geneChecker_error.log", "a"),
            )
            sort_gff.wait()

    # check genes (doublon ?)

    dico_genes = {}
    dico_gcount = {}
    c = len(glob.glob(pathOfFinalResults + "/*_genes_NT.fasta"))
    double_gene = 0
    if c == 1:
        for f in glob.glob(pathOfFinalResults + "/*_genes_NT.fasta"):
            f = f.split("_NT.fasta")[0]
            shutil.copy(
                f + "_AA.fasta",
                pathOfFinalResults + "/" + args.processName + "_final_genes_AA.fasta",
            )
            shutil.copy(
                f + "_NT.fasta",
                pathOfFinalResults + "/" + args.processName + "_final_genes_NT.fasta",
            )
    if c > 1 and args.numt == 0:
        for f in sorted(
            glob.glob(pathOfFinalResults + "/*_genes_NT.fasta"), key=os.path.getmtime
        ):
            for name, seq in read_fasta(open(f, "r")):
                gene = name.split("@")[1]
                if dico_genes.has_key(gene):
                    dico_gcount[gene] = dico_gcount.get(gene) + 1
                    if len(dico_genes.get(gene)) < len(seq):
                        dico_genes[gene] = seq
                else:
                    dico_genes[gene] = seq
                    dico_gcount[gene] = 1

        for key, value in dico_gcount.items():
            if value > 1 and c > 1:
                print(
                    "WARNING: "
                    + key
                    + " has been found more than once ("
                    + str(value)
                    + ") in the different mitochondrial contigs."
                )
                print(
                    "MitoFinder selected the longest sequence as the final sequence.\n"
                )
                logfile.write(
                    "WARNING : "
                    + key
                    + " has been found more than once ("
                    + str(value)
                    + ") in the different mitochondrial contigs."
                    + "\n"
                    + "Mitofinder selected the longest sequence as the final sequence.\n"
                    + "\n"
                )
                double_gene += 1

        final_fasta = open(
            pathOfFinalResults + "/" + args.processName + "_final_genes_NT.fasta", "w"
        )

        for key, value in dico_genes.items():
            final_fasta.write(">" + args.processName + "@" + key + "\n" + value + "\n")
        final_fasta.close()

        dico_genes = {}
        for f in sorted(
            glob.glob(pathOfFinalResults + "/*_genes_AA.fasta"), key=os.path.getmtime
        ):
            for name, seq in read_fasta(open(f, "r")):
                gene = name.split("@")[1]
                if dico_genes.has_key(gene):
                    if len(dico_genes.get(gene)) < len(seq):
                        dico_genes[gene] = seq
                else:
                    dico_genes[gene] = seq

        final_fasta = open(
            pathOfFinalResults + "/" + args.processName + "_final_genes_AA.fasta", "w"
        )

        for key, value in dico_genes.items():
            final_fasta.write(">" + args.processName + "@" + key + "\n" + value + "\n")
        final_fasta.close()

    print("## Final sequence saved to " + pathOfFinalResults)

    logfile.write("## Final sequence saved to " + pathOfFinalResults + "\n" + "\n")

    if double_gene == 1:
        print(
            "\n/!\\ WARNING /!\\ "
            + str(double_gene)
            + ' gene was found more than once suggesting either fragmentation, NUMT annotations, or potential contamination of your sequencing data.\nDifferent contigs may be part of different organisms thus "'
            + args.processName
            + '_final_genes_NT.fasta"'
            + ' and "'
            + args.processName
            + '_final_genes_AA.fasta" could be erroneous.\nWe recommend to check contigs and associated genes separately.\n'
        )
        logfile.write(
            "\n/!\\ WARNING /!\\ "
            + str(double_gene)
            + ' gene was found more than once suggesting either fragmentation, NUMT annotations, or potential contamination of your sequencing data.\nDifferent contigs may be part of different organisms thus "'
            + args.processName
            + '_final_genes_NT.fasta"'
            + ' and "'
            + args.processName
            + '_final_genes_AA.fasta" could be erroneous.\nWe recommend to check contigs and associated genes separately.\n'
        )
    elif double_gene > 1:
        print(
            "\n/!\\ WARNING /!\\ "
            + str(double_gene)
            + ' genes were found more than once suggesting either fragmentation, NUMT annotations, or potential contamination of your sequencing data.\nDifferent contigs may be part of different organisms thus "'
            + args.processName
            + '_final_genes_NT.fasta"'
            + ' and "'
            + args.processName
            + '_final_genes_AA.fasta" could be erroneous.\nWe recommend to check contigs and associated genes separately.\n'
        )
        logfile.write(
            "\n/!\\ WARNING /!\\ "
            + str(double_gene)
            + ' genes were found more than once suggesting either fragmentation, NUMT annotations, or potential contamination of your sequencing data.\nDifferent contigs may be part of different organisms thus "'
            + args.processName
            + '_final_genes_NT.fasta"'
            + ' and "'
            + args.processName
            + '_final_genes_AA.fasta" could be erroneous.\nWe recommend to check contigs and associated genes separately.\n'
        )

    if args.gap == 1 and args.numt == 0:
        print(
            '\n/!\\ WARNING /!\\ You have chosen to use the "allow-intron" option. Given that the annotation of exons is only based on similarity (BLASTX), we highly recommend to double-check the annotation.\n'
        )
        logfile.write(
            '\n/!\\ WARNING /!\\ You have chosen to use the "allow-intron" option. Given that the annotation of exons is only based on similarity (BLASTX), we highly recommend to double-check the annotation.\n'
        )

    if args.gap == 1 and args.numt == 1:
        print(
            "\n/!\\ WARNING /!\\ You have chosen to search for both introns and numts sequences at the same. This is difficult task. In that sense, please be cautious about the results! Some options are limited, for example, nWalk = 0 so you need to have a pretty good reference to optimize the annotation step."
        )
        logfile.write(
            "\n/!\\ WARNING /!\\ You have chosen to search for both introns and numts sequences at the same. This is difficult task. In that sense, please be cautious about the results! Some options are limited, for example, nWalk = 0 so you need to have a pretty good reference to optimize the annotation step."
        )

    # Cleaning

    trna_folder = pathtowork + "/" + args.processName + "_" + tRNA
    if os.path.exists(pathtowork + "/" + args.processName + "_" + tRNA):
        shutil.rmtree(pathtowork + "/" + args.processName + "_" + tRNA)
    if not os.path.exists(pathtowork + "/" + args.processName + "_" + tRNA):
        os.makedirs(pathtowork + "/" + args.processName + "_" + tRNA)
    if tRNA == "arwen":
        for f in glob.glob(pathOfFinalResults + "*.arwen"):
            shutil.copy(f, trna_folder + "/")
            os.remove(f)
        shutil.copy(pathOfFinalResults + "ARWEN.log", trna_folder + "/")
        os.remove(pathOfFinalResults + "ARWEN.log")
    elif tRNA == "trnascan":
        for f in glob.glob(pathOfFinalResults + "*.trnascan"):
            shutil.copy(f, trna_folder + "/")
            os.remove(f)
        shutil.copy(pathOfFinalResults + "tRNAscan-SE.log", trna_folder + "/")
        os.remove(pathOfFinalResults + "tRNAscan-SE.log")
    elif tRNA == "mitfi":
        for f in glob.glob(pathOfFinalResults + "*.mitfi"):
            shutil.copy(f, trna_folder + "/")
            os.remove(f)
        shutil.copy(pathOfFinalResults + "MiTFi.log", trna_folder + "/")
        os.remove(pathOfFinalResults + "MiTFi.log")

    tmpfiles = pathtowork + "/" + args.processName + "_tmp"
    if os.path.exists(tmpfiles):
        shutil.rmtree(tmpfiles)
    if not os.path.exists(tmpfiles):
        os.makedirs(tmpfiles)

    for f in glob.glob(pathtowork + "/ref*fasta.*"):
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + "/*.fasta.n*"):
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + "/*cmsearch*"):
        os.remove(f)
    for f in glob.glob(pathtowork + "/*.fasta.n*"):
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + "/important_fea*.fasta.*"):
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + "*_raw.gff"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*blast*"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*database*"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*.log"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/ref_for*"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*partial*"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*contig*.fasta"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathtowork + "/*.scafSeq.n*"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(tmpfiles + "/" + args.processName + "*.fasta"):
        os.remove(f)
    os.remove("genes_list")
    if args.genbk == False:
        for f in glob.glob(pathOfFinalResults + args.processName + "*.gb"):
            os.remove(f)
    for f in glob.glob(pathOfFinalResults + args.processName + "*.xml"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + args.processName + "*ref.fasta"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)
    for f in glob.glob(pathOfFinalResults + args.processName + "*ref.cds.fasta"):
        shutil.copy(f, tmpfiles + "/")
        os.remove(f)

    time = datetime.now() - start_time
    print("Total wall-clock time used by MitoFinder = " + str(time))
    logfile.write("\nTotal wall-clock time used by MitoFinder = " + str(time) + "\n")


if __name__ == "__main__":
    main()
