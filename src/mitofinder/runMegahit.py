import mitofinder.FirstBuildChecker as FirstBuildChecker
from subprocess import Popen
import subprocess
import time
import shlex
import os
import shutil
from shutil import copyfile

# deal with combo paired ans single reads
# fix path to app


def runMegahit(
    processName="teste",
    shortestContig=100,
    inputFile="teste.input",
    processorsToUse=4,
    refSeqFile=None,
    organismType=2,
    maxMemory="",
    logfile="logfile",
    override=False,
):
    pathToMegahit = megahitFolder
    bestBuild = None

    print("Starting Assembly step with MEGAHIT ")
    logfile = open(logfile, "a")
    logfile.write("Starting Assembly step with MEGAHIT " + "\n")

    pathToWork = os.getcwd() + "/"
    print("Result files will be saved here: ")
    print(pathToWork + processName + "_megahit/")
    logfile.write(
        "Result files will be saved here: "
        + "\n"
        + pathToWork
        + processName
        + "_megahit/"
        + "\n"
    )

    # copy input file to secondary folder to keep it organized
    if "/" in inputFile:
        destFile = pathToWork + str(inputFile.split("/")[-1])
    else:
        destFile = pathToWork + "/" + inputFile
    shutil.copyfile(inputFile, destFile)

    #####################################
    ######	    Run Megahit!!!     ######
    #####################################

    # create Megahit logfile:
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
    out = processName + "_megahit"
    megahit = "yes"
    if os.path.isdir(out) and override == False:
        print("\n####################################")
        print(
            "\n WARNING : "
            + pathToWork
            + out
            + " already exists. (use --override option)"
        )
        print("Mitofinder will skip MEGAHIT step")
        print(
            "\nIf you want to run MEGAHIT again, kill the mitofinder process, remove (or use --override) or rename the MEGAHIT result folder, and restart mitofinder\n"
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
            + "Mitofinder will skip MEGAHIT step"
            + "\n"
            + "\nIf you want to run MEGAHIT again, kill the mitofinder process, remove (or use --override) or rename the MEGAHIT result folder, and restart mitofinder\n"
            + "\n"
            + "#####################################\n"
            + "\n"
        )
        time.sleep(2)
        megahit = "no"
    elif os.path.isdir(out) and override == True:
        shutil.rmtree(out)
    if megahit == "yes":
        with open(pathToWork + "megahit.log", "w") as megahitLogFile:
            if t == "PE":
                if maxMemory == "":
                    command = (
                        "%smegahit -1 %s -2 %s -o %s --out-prefix %s --min-contig-len %s -t %s"
                        % (
                            pathToMegahit,
                            read1,
                            read2,
                            out,
                            out,
                            shortestContig,
                            processorsToUse,
                        )
                    )
                else:
                    command = (
                        "%smegahit -1 %s -2 %s -o %s --out-prefix %s --min-contig-len %s -t %s -m %s000000000"
                        % (
                            pathToMegahit,
                            read1,
                            read2,
                            out,
                            out,
                            shortestContig,
                            processorsToUse,
                            maxMemory,
                        )
                    )
                megahit = Popen(
                    command, stdout=megahitLogFile, stderr=megahitLogFile, shell=True
                )
                megahit.wait()
                if (
                    not os.path.isfile(
                        pathToWork + "/" + out + "/" + out + ".contigs.fa"
                    )
                    == True
                ):
                    print("\n ERROR: MEGAHIT didn't run well")
                    print("Please check log file : " + pathToWork + "megahit.log")
                    logfile.write(
                        "\n ERROR: MEGAHIT didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "megahit.log"
                        + "\n"
                    )
                    exit()
                # copyfile(pathToWork+"/"+out+"/"+out+".contigs.fa", pathToWork+"/"+processName+".scafSeq")
                # check Megahit output to see if reference sequence was built

            if t == "SE":
                if maxMemory == "":
                    command = (
                        "%smegahit -r %s -o %s --out-prefix %s --min-contig-len %s -t %s"
                        % (
                            pathToMegahit,
                            read1,
                            processName + "_megahit",
                            processName + "_megahit",
                            shortestContig,
                            processorsToUse,
                        )
                    )
                else:
                    command = (
                        "%smegahit -r %s -o %s --out-prefix %s --min-contig-len %s -t %s -m %s000000000"
                        % (
                            pathToMegahit,
                            read1,
                            processName + "_megahit",
                            processName + "_megahit",
                            shortestContig,
                            processorsToUse,
                            maxMemory,
                        )
                    )
                megahit = Popen(
                    command, stdout=megahitLogFile, stderr=megahitLogFile, shell=True
                )
                megahit.wait()
                if (
                    not os.path.isfile(
                        pathToWork + "/" + out + "/" + out + ".contigs.fa"
                    )
                    == True
                ):
                    print("\n MEGAHIT didn't run well")
                    print("Please check log file : " + pathToWork + "megahit.log")
                    logfile.write(
                        "\n ERROR: MEGAHIT didn't run well"
                        + "\n"
                        + "Please check log file : "
                        + pathToWork
                        + "megahit.log"
                        + "\n"
                    )
                    exit()
                # copyfile(pathToWork+"/"+out+"/"+out+".contigs.fa", pathToWork+"/"+processName+".scafSeq")

    logfile.close()
