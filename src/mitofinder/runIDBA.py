import mitofinder.FirstBuildChecker as FirstBuildChecker

from shutil import copyfile
from subprocess import Popen
import gzip
import os
import shlex
import shutil
import subprocess
import time
import logging

# Deal with combo single and paired reads
# fix app path


def runIDBA(
    processName="teste",
    shortestContig=100,
    inputFile="teste.input",
    processorsToUse=4,
    refSeqFile=None,
    organismType=2,
    override=False,
):
    bestBuild = None

    logging.info("Starting Assembly step with IDBA-UD ")

    pathToWork = os.getcwd() + "/"

    logging.info(
        "Result files will be saved here: "
        + "\n"
        + pathToWork
        + processName
        + "_idba/"
        + "\n"
    )

    # copy input file to secondary folder to keep it organized
    if "/" in inputFile:
        destFile = pathToWork + str(inputFile.split("/")[-1])
    else:
        destFile = pathToWork + "/" + inputFile
    shutil.copyfile(inputFile, destFile)

    #####################################
    ######	    Run IDBA-UD!!!     ######
    #####################################

    out = processName + "_idba"
    idba = "yes"
    if os.path.isdir(out) and override == False:
        logging.info(
            "\n####################################"
            + "\n"
            + "\n WARNING : "
            + pathToWork
            + out
            + " already exists. (use --override option)"
            + "\n"
            + "Mitofinder will skip idba step"
            + "\n"
            + "\nIf you want to run idba again, kill the mitofinder process, remove (or use --override) or rename the idba result folder, and restart mitofinder\n"
            + "\n"
            + "#####################################\n"
            + "\n"
        )
        time.sleep(2)
        idba = "no"
    elif os.path.isdir(out) and override == True:
        shutil.rmtree(out)
    # create IDBA logfile:
    if idba == "yes":
        with open(inputFile, "r") as InputFile:
            with open(pathToWork + "idba.log", "w") as idbaLogFile:
                for line in InputFile:
                    if "#" != line[0] and line != "\n":
                        configPart = (
                            line.lower()
                            .replace("\n", "")
                            .replace(" ", "")
                            .split("=")[0]
                        )
                        if configPart == "type":
                            t = line.replace("\n", "").replace(" ", "").split("=")[-1]
                        if configPart == "q1":
                            read1 = (
                                line.replace("\n", "").replace(" ", "").split("=")[-1]
                            )
                            if not "/" in read1:
                                read1 = "../" + read1
                            if read1[-3:] == ".gz":
                                with gzip.open(read1, "rb") as f_in:
                                    with open(read1[0:-3], "wb") as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                read1 = read1[0:-3]
                        if configPart == "q2":
                            read2 = (
                                line.replace("\n", "").replace(" ", "").split("=")[-1]
                            )
                            if not "/" in read2:
                                read2 = "../" + read2
                            if read2[-3:] == ".gz":
                                with gzip.open(read2, "rb") as f_in:
                                    with open(read2[0:-3], "wb") as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                read2 = read2[0:-3]

        with open(pathToWork + "idba.log", "a") as idbaLogFile:
            if t == "PE":
                print("Paired-end")
                logging.info("Paired-end" + "\n")
                read = processName + "_idba_read.fasta"
                command = "fq2fa --merge --filter %s %s %s" % (
                    read1,
                    read2,
                    read,
                )
                logging.info("Preparing data for IDBA-UD assembly")
                fq2fa = Popen(
                    command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True
                )
                fq2fa.wait()
                if not os.path.isfile(pathToWork + "/" + read) == True:
                    logging.error(
                        "\n ERROR: IDBA-UD didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "idba.log"
                        + "\n"
                    )
                    exit(1)
                command = "idba -r %s -o %s --num_threads %s" % (
                    read,
                    out,
                    processorsToUse,
                )
                print("Running assembly")
                idba = Popen(
                    command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True
                )
                idba.wait()
                # copyfile(pathToWork+"/"+out+"/contig.fa", pathToWork+"/"+processName+".scafSeq")

            if t == "SE":
                logging.info("Single-end")
                read = processName + "_idba_read.fasta"
                command = "fq2fa --filter %s %s" % (read1, read)
                logging.info("Preparing data for IDBA-UD assembly")
                fq2fa = Popen(
                    command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True
                )
                fq2fa.wait()
                if not os.path.isfile(pathToWork + "/" + read) == True:
                    logging.error(
                        "\n ERROR: IDBA-UD didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "idba.log"
                        + "\n"
                    )
                    exit(1)
                command = "idba -r %s -o %s --num_threads %s" % (
                    read,
                    out,
                    processorsToUse,
                )
                logging.info("Running assembly")
                idba = Popen(
                    command, stdout=idbaLogFile, stderr=idbaLogFile, shell=True
                )
                idba.wait()
                # copyfile(pathToWork+"/"+out+"/contig.fa", pathToWork+"/"+processName+".scafSeq")

        with gzip.open(read + ".gz", "wb") as f:
            f.write(read)
            os.remove(read)
