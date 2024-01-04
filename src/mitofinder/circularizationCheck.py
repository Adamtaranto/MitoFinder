from mitofinder.utils import check_files

from Bio import SeqIO, SearchIO

from subprocess import Popen
import shlex, sys, os
import logging


def cleanup_blast(blastTemp):
    if os.path.exists(blastTemp):
        logging.info(f"Cleanup blast output: {blastTemp}")
        os.remove(blastTemp)


def circularizationCheck(resultFile, circularSize, circularOffSet, tempdir):
    """
    Check, with blast, if there is a match between the start and the end of a sequence.
    Returns a tuple with (True, start, end) or False, accordingly.
    """
    check_files([resultFile])

    if not tempdir:
        tempdir = os.getcwd()
    elif not os.path.exists(tempdir):
        logging.warning(f"tmp dir {tempdir} not found. Write to cwd instead.")
        tempdir = os.getcwd()
    else:
        logging.info(f"Writing temp blast files to: {tempdir}")

    refSeq = SeqIO.read(resultFile, "fasta")
    sizeOfSeq = len(refSeq)

    try:
        command = "makeblastdb -in " + resultFile + " -dbtype nucl"
        args = shlex.split(command)
        logging.info(f"Building blast db: {command}")
        formatDB = Popen(args, stdout=open(os.devnull, "wb"))
        formatDB.wait()
        logging.info("makeblastdb complete.")
    except:
        logging.warning("formatDB during circularization check failed.")

        return (False, -1, -1)

    blastoutpath = os.path.join(tempdir, "circularization_check.blast.xml")
    with open(blastoutpath, "w") as blastResultFile:
        command = (
            "blastn -task blastn -db "
            + resultFile
            + " -query "
            + resultFile
            + " -outfmt 5"
        )  # call BLAST with XML output
        args = shlex.split(command)
        logging.info(f"Run BLAST search: {command}")
        blastAll = Popen(args, stdout=blastResultFile)
        blastAll.wait()
        logging.info("BLAST complete.")

    # Load hits
    blastparse = SearchIO.parse(blastoutpath, "blast-xml")

    """
	Let's loop through all blast results and see if there is a circularization.
	Do it by looking at all HSPs in the parse and see if there is an alignment of the ending of the sequence 
    with the start of that same sequence. It should have a considerable size, you don't want to say it circularized
	if only a couple of bases matched.
	Returns True or False, x_coordinate, y_coordinate
	x coordinate = starting point of circularization match
	y coordinate = ending point of circularization match
	"""
    for qresult in blastparse:  # in each query...
        for (
            hsp
        ) in (
            qresult.hsps
        ):  # loop through all HSPs looking for a circularization (perceived as a hsp with start somewhat close to the query finish)
            if (
                (hsp.query_range[0] >= 0 and hsp.query_range[0] <= circularOffSet)
                and (
                    hsp.hit_range[0] >= sizeOfSeq - hsp.aln_span - circularOffSet
                    and hsp.hit_range[0] <= sizeOfSeq + circularOffSet
                )
                and hsp.aln_span >= circularSize
                and hsp.aln_span < sizeOfSeq * 0.90
            ):
                if hsp.hit_range[0] < hsp.query_range[0]:
                    cleanup_blast(blastoutpath)
                    return (
                        True,
                        hsp.hit_range[0],
                        hsp.hit_range[1],
                    )  # it seems to have circularized, return True
                else:
                    cleanup_blast(blastoutpath)
                    return (True, hsp.query_range[0], hsp.query_range[1])

    # no circularization was observed in the for loop, so we exited it, just return false
    cleanup_blast(blastoutpath)
    return (False, -1, -1)


def main():
    # Config logging which called as main.
    logging.basicConfig(
        level=0, format="%(asctime)s:%(levelname)s:%(module)s:%(message)s"
    )

    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(
            "Usage: python -m mitofinder.circularizationCheck [fasta_file] [CircularSize] [CircularOffset] [tempdir] "
        )
    else:
        try:
            tempdir = sys.argv[4]
        except:
            tempdir = os.getcwd()
        print(
            circularizationCheck(
                sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), tempdir
            )
        )


if __name__ == "__main__":
    main()
