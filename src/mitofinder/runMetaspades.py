import mitofinder.FirstBuildChecker as FirstBuildChecker

from shutil import copyfile
from subprocess import Popen
import logging
import os
import shlex
import shutil
import subprocess
import time

# deal with combo paired and single reads


def runMetaspades(
    processName="teste",
    shortestContig=100,
    inputFile="teste.input",
    processorsToUse=4,
    refSeqFile=None,
    organismType=2,
    maxMemory="",
    override=False,
):
    bestBuild = None
    logging.info("Starting Assembly step with MetaSPAdes ")

    pathToWork = os.getcwd() + "/"
    logging.info(
        "Result files will be saved here: "
        + "\n"
        + pathToWork
        + processName
        + "_metaspades/"
        + "\n"
    )
    # copy input file to secondary folder to keep it organized
    if "/" in inputFile:
        destFile = pathToWork + str(inputFile.split("/")[-1])
    else:
        destFile = pathToWork + "/" + inputFile
    shutil.copyfile(inputFile, destFile)

    #####################################
    ######	 Run MetaSPAdes!!!     ######
    #####################################

    # create Metaspades logfile:

    with open(inputFile, "r") as InputFile:
        for line in InputFile:
            if "#" != line[0] and line != "\n":
                configPart = (
                    line.lower().replace("\n", "").replace(" ", "").split("=")[0]
                )
                if configPart == "type":
                    t = line.replace("\n", "").replace(" ", "").split("=")[-1]
                if configPart == "q1":
                    read1 = line.replace("\n", "").replace(" ", "").split("=")[-1]
                    if not "/" in read1:
                        read1 = "../" + read1
                if configPart == "q2":
                    read2 = line.replace("\n", "").replace(" ", "").split("=")[-1]
                    if not "/" in read2:
                        read2 = "../" + read2

    # try:
    out = processName + "_metaspades"
    metaspades = "yes"
    if os.path.isdir(out) and override == False:
        logging.warning(
            "####################################"
            + "\n WARNING : "
            + pathToWork
            + out
            + " already exists. (use --override option)"
            + "\n"
            + "Mitofinder will skip MetaSPAdes step"
            + "\n"
            + "\nIf you want to run MetaSPAdes again, kill the mitofinder process, remove (or use --override) or rename the MetaSPAdes result folder, and restart mitofinder\n"
            + "#####################################\n"
        )
        time.sleep(2)
        metaspades = "no"
    elif os.path.isdir(out) and override == True:
        shutil.rmtree(out)
    if metaspades == "yes":
        with open(pathToWork + "metaspades.log", "w") as metaspadesLogFile:
            if t == "PE":
                if maxMemory == "":
                    command = "metaspades.py -1 %s -2 %s -o %s -t %s" % (
                        read1,
                        read2,
                        out,
                        processorsToUse,
                    )
                else:
                    command = "metaspades.py -1 %s -2 %s -o %s -t %s -m %s" % (
                        read1,
                        read2,
                        out,
                        processorsToUse,
                        maxMemory,
                    )
                metaspades = Popen(
                    command,
                    stdout=metaspadesLogFile,
                    stderr=metaspadesLogFile,
                    shell=True,
                )
                metaspades.wait()
                if (
                    not os.path.isfile(pathToWork + "/" + out + "/" + "scaffolds.fasta")
                    == True
                ):
                    logging.error(
                        "\n ERROR: MetaSPAdes didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "metaspades.log"
                        + "\n"
                    )
                    exit(1)
                # copyfile(pathToWork+"/"+out+"/"+"scaffolds.fasta", pathToWork+"/"+processName+".scafSeq")
                # check Megahit output to see if reference sequence was built

            if t == "SE":
                if maxMemory == "":
                    command = "metaspades.py -s %s -o %s -t %s" % (
                        read1,
                        out,
                        processorsToUse,
                    )
                else:
                    command = "metaspades.py -s %s -o %s -t %s -m %s" % (
                        read1,
                        out,
                        processorsToUse,
                        maxMemory,
                    )
                metaspades = Popen(
                    command,
                    stdout=metaspadesLogFile,
                    stderr=metaspadesLogFile,
                    shell=True,
                )
                metaspades.wait()
                if (
                    not os.path.isfile(pathToWork + "/" + out + "/" + "scaffolds.fasta")
                    == True
                ):
                    logging.error(
                        "\n MetaSPAdes didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "metaspades.log"
                        + "\n"
                    )
                    exit(1)
                # copyfile(pathToWork+"/"+out+"/"+"scaffolds.fasta", pathToWork+"/"+processName+".scafSeq")
