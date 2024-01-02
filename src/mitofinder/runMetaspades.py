import mitofinder.FirstBuildChecker as FirstBuildChecker
from shutil import copyfile
from subprocess import Popen
import os
import shlex
import shutil
import subprocess
import time


def runMetaspades(
    processName="teste",
    shortestContig=100,
    inputFile="teste.input",
    processorsToUse=4,
    metaspadesFolder="installed",
    refSeqFile=None,
    organismType=2,
    blastFolder="installed",
    maxMemory="",
    logfile="logfile",
    override=False,
):
    pathToMetaspades = metaspadesFolder
    bestBuild = None

    logfile = open(logfile, "a")
    print("Starting Assembly step with MetaSPAdes ")
    logfile.write("Starting Assembly step with MetaSPAdes " + "\n")

    pathToWork = os.getcwd() + "/"
    print("Result files will be saved here: ")
    print(pathToWork + processName + "_metaspades/")
    logfile.write(
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
        print("\n####################################")
        print(
            "\n WARNING : "
            + pathToWork
            + out
            + " already exists. (use --override option)"
        )
        print("Mitofinder will skip MetaSPAdes step")
        print(
            "\nIf you want to run MetaSPAdes again, kill the mitofinder process, remove (or use --override) or rename the MetaSPAdes result folder, and restart mitofinder\n"
        )
        print("#####################################\n")
        logfile.write(
            "\n####################################"
            + "\n"
            + "\n WARNING : "
            + pathToWork
            + out
            + " already exists. (use --override option)"
            + "\n"
            + "Mitofinder will skip MetaSPAdes step"
            + "\n"
            + "\nIf you want to run MetaSPAdes again, kill the mitofinder process, remove (or use --override) or rename the MetaSPAdes result folder, and restart mitofinder\n"
            + "\n"
            + "#####################################\n"
            + "\n"
        )
        time.sleep(2)
        metaspades = "no"
    elif os.path.isdir(out) and override == True:
        shutil.rmtree(out)
    if metaspades == "yes":
        with open(pathToWork + "metaspades.log", "w") as metaspadesLogFile:
            if t == "PE":
                if maxMemory == "":
                    command = "%smetaspades.py -1 %s -2 %s -o %s -t %s" % (
                        pathToMetaspades,
                        read1,
                        read2,
                        out,
                        processorsToUse,
                    )
                else:
                    command = "%smetaspades.py -1 %s -2 %s -o %s -t %s -m %s" % (
                        pathToMetaspades,
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
                    print("\n ERROR: MetaSPAdes didn't run well")
                    print("Please check log file : " + pathToWork + "metaspades.log")
                    logfile.write(
                        "\n ERROR: MetaSPAdes didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "metaspades.log"
                        + "\n"
                    )
                    exit()
                # copyfile(pathToWork+"/"+out+"/"+"scaffolds.fasta", pathToWork+"/"+processName+".scafSeq")
                # check Megahit output to see if reference sequence was built

            if t == "SE":
                if maxMemory == "":
                    command = "%smetaspades.py -s %s -o %s -t %s" % (
                        pathToMetaspades,
                        read1,
                        out,
                        processorsToUse,
                    )
                else:
                    command = "%smetaspades.py -s %s -o %s -t %s -m %s" % (
                        pathToMetaspades,
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
                    print("\n MetaSPAdes didn't run well")
                    print("Please check log file : " + pathToWork + "metaspades.log")
                    logfile.write(
                        "\n MetaSPAdes didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "metaspades.log"
                        + "\n"
                    )
                    exit()
                # copyfile(pathToWork+"/"+out+"/"+"scaffolds.fasta", pathToWork+"/"+processName+".scafSeq")

    logfile.close()
