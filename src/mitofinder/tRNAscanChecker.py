from mitofinder.utils import is_avail, add_to_path

from Bio import SeqIO
from Bio.Data import CodonTable

from subprocess import Popen
import shlex, sys, os
import logging
import shutil


class Assembly:
    """
    Class to hold the sequence from best_query and it's tRNAscan-SE results, as well as features.
    """

    def __init__(
        self, resultFile, tRNAscanResultFile, hasCircularized, tRNAscan, organismType
    ):
        if resultFile[-6:] != ".fasta":
            self.refSeq = SeqIO.read(resultFile, "genbank")
        else:
            self.refSeq = SeqIO.read(resultFile, "fasta")
        self.tRNAs = []
        self.checkFeatures = None
        self.validContigs = 1
        checkedtRNAs = []
        self.hasCircularized = hasCircularized
        startCheck = False
        if tRNAscanResultFile != None:
            with open(tRNAscanResultFile, "r") as tRNAscanFile:
                if tRNAscanResultFile[-6:] == ".mitfi":
                    module_dir = os.path.dirname(__file__)
                    module_dir = os.path.abspath(os.path.join(module_dir, "../../bin"))
                    module_dir = os.path.join(module_dir, "mitfi")
                    dico_tRNA = {}
                    for line in open(
                        module_dir + "/Codon-Amino_Acid_Abbreviations.csv"
                    ):
                        dico_tRNA[line.rstrip().split(";")[3]] = line.split(";")[2]
                    c = 1
                    for line in tRNAscanFile:
                        if line[:7] == "#header":
                            """
                            look for the start of tRNA results as in:
                            #header	start	stop	evalue	AC	AA	model	strand
                            """
                            startCheck = True
                        elif (
                            startCheck == True and line[0] != "#"
                        ):  # we found the start, so let's start saving this info!
                            tRNAinfos = line.split()
                            # sequence name is tRNAinfos[0], but isn't really needed
                            tRNAnumber = c
                            c += 1
                            tRNAstrand = tRNAinfos[8]
                            if tRNAstrand == "+":
                                tRNAcoordinates = (int(tRNAinfos[1]), int(tRNAinfos[2]))
                            elif tRNAstrand == "-":
                                tRNAcoordinates = (int(tRNAinfos[2]), int(tRNAinfos[1]))
                            tRNAtype = str(tRNAinfos[6])
                            tRNAtype = dico_tRNA.get(tRNAtype[0])
                            if tRNAtype in checkedtRNAs:
                                tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
                            checkedtRNAs.append(tRNAtype)
                            tRNAcodon = tRNAinfos[5]
                            tRNAintronBegin = 0
                            tRNAintronEnd = 0
                            if tRNAintronBegin > 0:
                                logging.warning(
                                    f"WARNING: {tRNAtype} found with an intron. Strange result, check .trnascan file"
                                )
                            tRNAscore = float(tRNAinfos[3])
                            thisRNA = self.tRNA(
                                self,
                                tRNAnumber,
                                tRNAcoordinates,
                                tRNAtype,
                                tRNAcodon,
                                tRNAscore,
                                tRNAintronBegin,
                                tRNAintronEnd,
                            )
                            self.tRNAs.append(thisRNA)
                elif tRNAscanResultFile[-9:] == ".trnascan":
                    for line in tRNAscanFile:
                        if line[0] == "-":
                            """
                            look for the start of tRNA results as in:
                            tRNA #	Begin	End  	Type	Codon	Begin	End	Score
                            """
                            startCheck = True
                        elif (
                            startCheck == True
                        ):  # we found the start, so let's start saving this info!
                            tRNAinfos = line.split()
                            # sequence name is tRNAinfos[0], but isn't really needed
                            tRNAnumber = int(tRNAinfos[1])
                            tRNAcoordinates = (int(tRNAinfos[2]), int(tRNAinfos[3]))
                            tRNAtype = tRNAinfos[4]
                            if tRNAtype in checkedtRNAs:
                                tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
                            checkedtRNAs.append(tRNAinfos[4])
                            tRNAcodon = tRNAinfos[5]
                            tRNAintronBegin = int(tRNAinfos[6])
                            tRNAintronEnd = int(tRNAinfos[7])
                            if tRNAintronBegin > 0:
                                logging.warning(
                                    "	WARNING: %s found with an intron. Strange result, check .trnascan file"
                                    % tRNAtype
                                )
                            tRNAscore = float(tRNAinfos[8])
                            thisRNA = self.tRNA(
                                self,
                                tRNAnumber,
                                tRNAcoordinates,
                                tRNAtype,
                                tRNAcodon,
                                tRNAscore,
                                tRNAintronBegin,
                                tRNAintronEnd,
                            )
                            self.tRNAs.append(thisRNA)
                elif tRNAscanResultFile[-6:] == ".arwen":
                    for line in tRNAscanFile:
                        if "genes found" in line:
                            """
                            look for the start of tRNA results as in:
                            tRNA #	Begin	End  	Type	Codon	Begin	End	Score

                            arwen
                            # mtRNA-Type [Begin,End] N (codon)


                            """
                            startCheck = True
                        elif (
                            startCheck == True
                        ):  # we found the start, so let's start saving this info!
                            tRNAinfos = line.split()
                            # sequence name is tRNAinfos[0], but isn't really needed
                            tRNAnumber = int(line.split(" ")[0])
                            if "c[" in line:
                                tRNAcoordinates = (
                                    line.split("]")[0].split(",")[1],
                                    line.split("[")[1].split(",")[0],
                                )
                            else:
                                tRNAcoordinates = (
                                    line.split("[")[1].split(",")[0],
                                    line.split("]")[0].split(",")[1],
                                )
                            tRNAtype = line.split("-")[1].split(" ")[0]
                            if tRNAtype in checkedtRNAs:
                                tRNAtype += str(checkedtRNAs.count(tRNAtype) + 1)
                            checkedtRNAs.append(tRNAtype)
                            tRNAcodon = line.split("(")[1].split(")")[0].upper()
                            tRNAintronBegin = 0
                            tRNAintronEnd = 0
                            tRNAscore = 50.0
                            thisRNA = self.tRNA(
                                self,
                                tRNAnumber,
                                tRNAcoordinates,
                                tRNAtype,
                                tRNAcodon,
                                tRNAscore,
                                tRNAintronBegin,
                                tRNAintronEnd,
                            )
                            self.tRNAs.append(thisRNA)

    def __len__(self):
        return len(self.tRNAs)

    def isCircular(self):
        return self.hasCircularized

    class tRNA:
        """
        Class to hold the tRNAscan-SE info for each tRNA in Assembly()
        """

        def __init__(
            self,
            motherSeq,
            tRNAnumber,
            tRNAcoordinates,
            tRNAtype,
            tRNAcodon,
            tRNAscore,
            tRNAintronBegin,
            tRNAintronEnd,
        ):
            self.tRNAnumber = tRNAnumber
            reverseStart = False
            if int(tRNAcoordinates[0]) > int(tRNAcoordinates[1]):
                reverseStart = True
            tRNAstart = min(int(tRNAcoordinates[0]), int(tRNAcoordinates[1]))
            tRNAend = max(int(tRNAcoordinates[0]), int(tRNAcoordinates[1]))
            tRNAstart -= 1
            tRNAstart = max(0, tRNAstart)
            if reverseStart is False:
                tRNAcoordinates = [tRNAstart, tRNAend]
            else:
                tRNAcoordinates = [tRNAend, tRNAstart]
            self.tRNAcoordinates = tRNAcoordinates
            self.tRNAtype = tRNAtype
            self.tRNAcodon = tRNAcodon
            self.tRNAscore = tRNAscore
            self.tRNAintronBegin = tRNAintronBegin
            self.tRNAintronEnd = tRNAintronEnd
            self.motherSeq = motherSeq

        def __len__(self):
            return self.coordinates()[1] - self.coordinates()[0]

        def coordinates(self):
            return self.tRNAcoordinates

        def number(self):
            return self.tRNAnumber

        def typeOfRna(self):
            return self.tRNAtype

        def codon(self):
            return self.tRNAcodon

        def score(self):
            return self.tRNAscore


def tRNAscanCheck(
    resultFile=None,
    hasCircularized=False,
    skipTRNA=False,
    organismType=2,
    coveCutOff=7,
    buildBacteria=False,
    buildArchea=False,
    tRNAscan="mitfi",
):
    """
    Use tRNAscan-SE to look for tRNAs and hold it's positions and scores in the tRNA Class.
    """
    if not skipTRNA:
        module_dir = os.path.dirname(__file__)
        module_dir = os.path.abspath(os.path.join(module_dir, "../../bin"))
        scanInput = resultFile

        if tRNAscan == "mitfi":
            # If not avail on PATH add bundled package to path.
            if not is_avail(["mitfi.jar"], kill=False):
                module_dir = os.path.join(module_dir, "mitfi")
                logging.info(f"Fallback to bundled Mitfi: {module_dir}")
                add_to_path(module_dir)
                # or let mitofinder handle paths later
                MitFiFolder = module_dir
            else:
                MitFiFolder = os.path.dirname(shutil.which("mitfi.jar"))
            # Set output file extension
            outputName = resultFile[0:-6] + ".mitfi"

        elif tRNAscan == "trnascan":
            # Kill if not avail on PATH
            is_avail(["trnascan-SE"])
            # Set output file extension
            outputName = resultFile[0:-6] + ".trnascan"

        elif tRNAscan == "arwen":
            # If not avail on PATH add bundled package to path.
            if not is_avail(["arwen"], kill=False):
                module_dir = os.path.join(module_dir, "arwen")
                logging.info(f"Fallback to bundled Arwen: {module_dir}")
                add_to_path(module_dir)
                ArwenFolder = module_dir
            else:
                logging.info(f'Using arwen from: {shutil.which("arwen")}')
                ArwenFolder = os.path.dirname(shutil.which("arwen"))
                # or let mitofinder handle paths later
            # Set output file extension
            outputName = resultFile[0:-6] + ".arwen"

        # Adding the Arwen folder to these environments in the OS to avoid errors being thrown by Arwen
        if tRNAscan == "mitfi":
            try:
                os.environ["PATH"] += os.pathsep + MitFiFolder
            except KeyError:
                pass
            except:
                logging.error("Error setting MitFi PATH.")

            try:
                os.environ["PERL5LIB"] += os.pathsep + MitFiFolder
            except KeyError:
                pass
            except:
                logging.error("Error setting PERL5LIB MitFi path.")
        elif tRNAscan == "arwen":
            try:
                os.environ["PATH"] += os.pathsep + ArwenFolder
            except KeyError:
                pass
            except:
                logging.error("Error setting Arwen path.")

            try:
                os.environ["PERL5LIB"] += os.pathsep + ArwenFolder
            except KeyError:
                pass
            except:
                logging.error("Error setting Arwen PERL5LIB path.")

        if os.path.exists(
            outputName
        ):  # remove result file if it already exists so that tRNAscan doesn't throw another error
            os.remove(outputName)

        if tRNAscan == "mitfi":
            out = outputName
            try:
                with open("MiTFi.log", "w") as tRNAscanLog:
                    command = (
                        "java -jar "
                        + MitFiFolder
                        + "/mitfi.jar -code "
                        + str(organismType)
                        + " "
                        + scanInput
                    )
                    args = shlex.split(command)
                    logging.info(f"Calling: {" ".join(args)}")
                    tRNAscanRun = Popen(
                        args, stdout=open(outputName, "w"), stderr=tRNAscanLog
                    )  # cwd=MitFiFolder
                    tRNAscanRun.wait()

                thisSequenceResult = Assembly(
                    resultFile, outputName, hasCircularized, tRNAscan, organismType
                )
                return thisSequenceResult
            except:
                if os.path.exists(
                    out
                ):  # remove result file if it already exists so that tRNAscan doesn't throw another error
                    os.remove(out)
                logging.warning("MiTFi failed.")
                return False

        elif tRNAscan == "trnascan":
            # down here, we check organism type to make tRNAscan-SE run in appropriate mode
            if buildBacteria == True:
                organismFlag = "-B "  # bacterial
            elif buildArchea == True:
                organismFlag = "-A "  # archea
            elif organismType == 1:
                organismFlag = ""  # leave it empty to use eukariotyc search
            else:
                organismFlag = "-O "

            # use different genetic code?
            if organismType == 2:  # vertebrate
                geneticCode = "-g lib/gcode/gcode.vertmito "
            elif organismType == 3:  # yeast
                geneticCode = "-g lib/gcode/gcode.ystmito "
            elif organismType == 4:  # mold, protozoan
                geneticCode = "-g lib/gcode/gcode.othmito "
            elif organismType == 5:  # invertebrate
                geneticCode = "-g lib/gcode/gcode.invmito "
            elif organismType == 6:  # ciliate
                geneticCode = "-g lib/gcode/gcode.cilnuc "
            elif organismType == 9:  # echinodermata
                geneticCode = "-g lib/gcode/gcode.echdmito "
            else:
                geneticCode = ""

            try:
                with open("tRNAscan-SE.log", "w") as tRNAscanLog:
                    command = (
                        "tRNAscan-SE -X "
                        + str(coveCutOff)
                        + " "
                        + geneticCode
                        + organismFlag
                        + "-o "
                        + outputName
                        + " "
                        + scanInput
                    )
                    args = shlex.split(command)
                    logging.info(f"Calling: {command}")
                    tRNAscanRun = Popen(args, stdout=tRNAscanLog, stderr=tRNAscanLog)
                    tRNAscanRun.wait()

                thisSequenceResult = Assembly(
                    resultFile, outputName, hasCircularized, tRNAscan, organismType
                )
                return thisSequenceResult
            except:
                logging.warning("tRNAscan-SE failed. Check it's logs for more details.")
                return False

        elif tRNAscan == "arwen":
            logging.info("Preparing Arwen command.")
            # down here, we check organism type to make ARWEN run in appropriate mode
            if buildBacteria == True:
                geneticCode = "-gcbact "  # bacterial
            elif buildArchea == True:
                geneticCode = "-gcbact "  # archea
            elif organismType == 1:
                geneticCode = "-gcstd "  # standard
            elif organismType == 2:  # vertebrate
                geneticCode = "-gcvert "
            elif organismType == 3:  # yeast
                geneticCode = "-gcyeast "
            elif organismType == 4:  # mold, protozoan
                geneticCode = "-gcprot "
            elif organismType == 5:  # invertebrate
                geneticCode = "-gcinvert "
            elif organismType == 6:  # ciliate
                geneticCode = "-gcciliate "
            elif organismType == 9:  # echinodermata
                geneticCode = "-gcflatworm "
            else:
                geneticCode = ""

            try:
                with open("ARWEN.log", "w") as tRNAscanLog:
                    command = (
                        ArwenFolder
                        + "/arwen "
                        + " "
                        + geneticCode
                        + "-o "
                        + outputName
                        + " -w "
                        + scanInput
                    )
                    args = shlex.split(command)
                    logging.info(f"Calling: {command}")
                    tRNAscanRun = Popen(
                        args, stdout=tRNAscanLog, stderr=tRNAscanLog
                    )  # cwd=ArwenFolder
                    tRNAscanRun.wait()

                thisSequenceResult = Assembly(
                    resultFile, outputName, hasCircularized, tRNAscan, organismType
                )
                return thisSequenceResult
            except:
                logging.warning("Arwen failed. Check log files for more details.")
                return False
    else:  # if tRNA = True
        return Assembly(resultFile, None, hasCircularized)


def tRNAconvert(inputTRNA):
    """
    Gets as input: a list of tRNAs in tRNA-Phe format or TRNF
    Returns: a list or a string (if input was only a string of tRNA)
            with the format tRNA-Phe
    DOES NOT ALTER THE FINAL NAME, just the checks
    """
    returnString = False
    if type(inputTRNA) != list:
        inputTRNA = [str(inputTRNA)]
        returnString = True
    outputList = []
    for convertTRNA in inputTRNA:
        dictOftRNAs = {
            "trnf": "trna-phe",
            "trnv": "trna-val",
            "trnl": "trna-leu",
            "trni": "trna-ile",
            "trnq": "trna-gln",
            "trnm": "trna-met",
            "trnw": "trna-trp",
            "trna": "trna-ala",
            "trnn": "trna-asn",
            "trnc": "trna-cys",
            "trny": "trna-tyr",
            "trns": "trna-ser",
            "trnd": "trna-asp",
            "trnk": "trna-lys",
            "trng": "trna-gly",
            "trnr": "trna-arg",
            "trnh": "trna-his",
            "trne": "trna-glu",
            "trnt": "trna-thr",
            "trnp": "trna-pro",
        }

        thisTRNA = convertTRNA.lower()
        thisPartialTRNA = thisTRNA[0:4]
        if len(thisTRNA) <= 6:
            if thisPartialTRNA in dictOftRNAs:
                if len(thisTRNA) == 4:
                    finalName = dictOftRNAs[thisTRNA]
                else:
                    if thisTRNA[-1] == "1":
                        finalName = dictOftRNAs[thisPartialTRNA]
                    else:
                        finalName = dictOftRNAs[thisPartialTRNA] + thisTRNA[-1]
        else:
            finalName = thisTRNA

        outputList.append(finalName)

    if returnString == False:
        return outputList
    else:
        return outputList[0]


def prettyRNAName(tRNAName):
    """
    Just a cosmetic function to properly format tRNA results that come in lowercase to the
    tRNA-Phe format.
    """
    return (
        tRNAName.replace("trna", "tRNA")[0:5]
        + tRNAName[tRNAName.index("-") + 1 :].capitalize()
    )
