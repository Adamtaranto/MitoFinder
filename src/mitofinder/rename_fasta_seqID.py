import sys
import os.path
from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
import logging


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


def main():
    logging.basicConfig(
        level=0, format="%(asctime)s:%(levelname)s:%(module)s:%(message)s"
    )

    usemsg = "Usage: python -m rename_fasta_seqID [seqID] [mito_fasta] [outfile] [suffix (int)] [direction (+/-)] [rename (True/False)]"

    if len(sys.argv) == 1:
        print(usemsg)
        exit(1)
    elif sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(usemsg)
        exit(0)
    elif len(sys.argv) < 7:
        logging.warning(f"Found {len(sys.argv) - 1 } args. Expected 6.")
        print(usemsg)
        exit(1)

    # Parse args
    seqID = sys.argv[1]
    mito_fasta = sys.argv[2]
    outfile = str(sys.argv[3])
    suffix = str(sys.argv[4])
    direction = sys.argv[5]
    rename = sys.argv[6]

    if rename == "True":
        logging.info(f"Appending suffic .{suffix} to output sequence name.")
    else:
        logging.info("Sequence will not be renamed.")

    resultFile = SeqIO.read(open(mito_fasta, "r"), "fasta")

    logging.info(f"Writing processed sequence to: {outfile}")
    fout = open(outfile, "w")

    if direction == "+":
        if rename == "True":
            fout.write(">" + seqID + "." + suffix + "\n" + str(resultFile.seq) + "\n")
        else:
            fout.write(">" + str(resultFile.id) + "\n" + str(resultFile.seq) + "\n")
    elif direction == "-":
        logging.info("reverse complementing sequence")
        if rename == "True":
            fout.write(
                ">"
                + seqID
                + "."
                + suffix
                + "\n"
                + str(resultFile.seq.reverse_complement())
                + "\n"
            )
        else:
            fout.write(
                ">"
                + str(resultFile.id)
                + " (reverse)\n"
                + str(resultFile.seq.reverse_complement())
                + "\n"
            )


if __name__ == "__main__":
    main()
