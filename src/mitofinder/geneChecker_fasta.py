from mitofinder import genbankOutput, tRNAscanChecker
from mitofinder.tRNAscanChecker import tRNAconvert, prettyRNAName
from mitofinder.utils import is_avail

from Bio import SeqIO, SearchIO, SeqFeature

from subprocess import Popen
import os
import shlex
import shutil
import sys
import logging


class Alignment:
    """
    Class to hold the alignments from Blast+.
    Initially it was for Needle, so it was kept for backwards compatibility.
    """

    def __init__(self, seq1, seq2, alignLength):  # seqId1, seqId2
        self.seq1 = seq1
        self.seq2 = seq2
        self.alignLength = alignLength
        self.alignseq = ""
        self.startBase = 0
        self.endBase = 0
        self.frame = 1
        self.translationTable = 2
        self.refSeq = None
        self.seqFound = ""

    def __len__(self):
        return self.alignLength

    def __str__(self):
        return self.seqFound

    def __lt__(self, other):
        return self.startBase < other.startBase


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


def geneCheck(
    fastaReference=None,
    resultFile=None,
    cutoffEquality_prot=None,
    cutoffEquality_nucl=None,
    usedOwnGenBankReference=None,
    organismType=2,
    alignCutOff=45,
    out_blast=None,
    blasteVal=None,
):
    """
    Returns a tuple with 2 dictionaries, one with the features found and another with features to look for.
    """
    # 	record = SeqIO.read(genBankReference, "genbank")
    refSeq = SeqIO.read(resultFile, "fasta")
    listOfImportantFeatures = {}
    logging.info(
        "Checking protein-coding genes, tRNAs and rRNAs from reference with organismType=%s..."
        % organismType
    )
    listOfPresentFeatures = {}
    listOfSplits = []
    listOfCompleteGenes = []
    """
	Do protein-coding genes first!
	"""
    # let's create the fasta file and the list of features we are looking for
    refGenes = open(fastaReference)
    importantFeaturesFile = open("important_features.fasta", "w")
    nbrgene = 0
    for name, seq in read_fasta(refGenes):
        if name.split("@")[1] != "rrnL" and name.split("@")[1] != "rrnS":
            nbrgene += 1
            featureName = name.split("@")[1]
            importantFeaturesFile.write(">" + featureName + "\n" + seq + "\n")
            listOfImportantFeatures[featureName] = seq
    importantFeaturesFile.close()
    # 	with open('important_features.fasta', 'w') as importantFeaturesFile:
    # 		for feature in record.features:
    # 			if feature.type.lower() == 'cds':
    # 				if 'gene' in feature.qualifiers:
    # 					featureName = feature.qualifiers['gene'][0]
    # 				elif 'product' in feature.qualifiers:
    # 					featureName = feature.qualifiers['product'][0]
    # 				featureName = ''.join(featureName.split())
    # 				if featureName in listOfImportantFeatures:
    # 					featureName += '_' + str(listOfImportantFeatures.keys().count(featureName) + 1)
    #
    # 				importantFeaturesFile.write('>' + featureName + '\n')
    # 				if 'translation' in feature.qualifiers:
    # 					importantFeaturesFile.write(str(feature.qualifiers['translation'][0]) + '\n')
    # 				else:
    # 					importantFeaturesFile.write(str(feature.extract(record).seq.translate(table=organismType,to_stop=True))+'\n')
    # 					logging.info('		WARNING: reference did not give a CDS translation for %s. Creating our own from refSeq.' \)
    # 						% featureName
    # 				listOfImportantFeatures[featureName] = feature

    # running blast
    if nbrgene > 0:
        logging.info("Formatting database for blast...")
        command = "makeblastdb -in important_features.fasta -dbtype prot"  # need to formatdb refseq first

        args = shlex.split(command)
        formatDB = Popen(args, stdout=open(os.devnull, "wb"))
        formatDB.wait()

        # logging.info("Running blast against refSeq to determine if a hit was built...")
        with open("important_features.blast.xml", "w") as blastResultFile:
            if usedOwnGenBankReference == True:  # using a personal genbank reference
                command = (
                    "blastx -db important_features.fasta -query "
                    + resultFile
                    + " -evalue "
                    + str(blasteVal)
                    + " -outfmt 5 -num_threads 2 -query_gencode "
                    + str(organismType)
                    + " -seg no"
                )  # call BLAST with XML output
            else:  # using a non personal genbank reference
                logging.info("Genetic code: ", str(organismType))
                command = (
                    "blastx -db important_features.fasta -query "
                    + resultFile
                    + "-evalue "
                    + str(blasteVal)
                    + " -outfmt 5 -num_threads 2 -query_gencode "
                    + str(organismType)
                    + " -seg no"
                )  # call BLAST with XML output
            args = shlex.split(command)
            blastAll = Popen(args, stdout=blastResultFile)
            blastAll.wait()
        # SearchIO object handler and checker for best hit separation
        listOfSplits = []
        listOfCompleteGenes = []
        blastparse = SearchIO.parse(
            "important_features.blast.xml", "blast-xml"
        )  # get all queries
        listOfPresentFeatures = {}
        for qresult in blastparse:  # in each query, let's look for a good hit
            for qhit in qresult.hits:
                for (
                    hsp
                ) in qhit.hsps:  # hsp object checking, this contains the alignment info
                    featureName = qhit.id
                    if float(str(hsp.ident_num) + ".00") / float(
                        str(hsp.aln_span) + ".00"
                    ) * 100 >= float(cutoffEquality_prot):
                        if featureName in listOfImportantFeatures:
                            targetFeature = listOfImportantFeatures[featureName]
                            if (
                                hsp.aln_span * 3
                                >= (len(targetFeature * 3) + 3) * alignCutOff / 100
                            ):
                                startBase = (
                                    min(hsp.query_range[0], hsp.query_range[1]) + 1
                                )
                                endBase = max(hsp.query_range[0], hsp.query_range[1])
                                alignLen = (endBase - startBase) + 1
                                if alignLen >= (len(targetFeature * 3) + 3) * 0.10:
                                    featureFrame = hsp.query_frame
                                    seqName = featureName
                                    alignment = Alignment(
                                        featureName, seqName, alignLen
                                    )
                                    alignment.refSeq = refSeq
                                    alignment.translationTable = organismType
                                    alignment.frame = featureFrame
                                    alignment.startBase = startBase
                                    alignment.endBase = endBase
                                    alignment.seqFound = refSeq.seq[
                                        (startBase - 1) : endBase
                                    ]
                                    listOfPresentFeatures[featureName] = (
                                        listOfImportantFeatures[qhit.id],
                                        alignment,
                                        featureFrame <= -1,
                                    )
                                    if alignLen >= (len(targetFeature * 3) + 3) * 0.99:
                                        # if we've already built a lot, dont even bother with finding splits
                                        listOfCompleteGenes.append(featureName)
                                    break

        # exit()
        # copying the blast result in order for this info to be assessed later if the user desires
        shutil.copyfile(
            "important_features.blast.xml", out_blast + "_ref.cds.blast.xml"
        )
        os.remove("important_features.blast.xml")
        shutil.copyfile("important_features.fasta", out_blast + "_ref.cds.fasta")
    os.remove("important_features.fasta")

    # let's create the fasta file and the list of features we are looking for
    refGenes = open(fastaReference)
    importantFeaturesFile = open("important_features.fasta", "w")
    nbrRNA = 0
    for name, seq in read_fasta(refGenes):
        if name.split("@")[1] == "rrnL" or name.split("@")[1] == "rrnS":
            nbrRNA += 1
            featureName = name.split("@")[1]
            importantFeaturesFile.write(">" + featureName + "\n" + seq + "\n")
            listOfImportantFeatures[featureName] = seq
    importantFeaturesFile.close()

    # 	with open('important_features.fasta', 'w') as importantFeaturesFile:
    # 		for feature in record.features:
    # 			if feature.type == 'rRNA' :
    # 				if 'gene' in feature.qualifiers:
    # 					featureName = feature.qualifiers['gene'][0]
    # 					featureName = ''.join(featureName.split())
    # 				elif 'product' in feature.qualifiers:
    # 					featureName = feature.qualifiers['product'][0]
    # 					featureName = ''.join(featureName.split())
    # 				if featureName in listOfImportantFeatures:
    # 					featureName += str(listOfImportantFeatures.keys().count(featureName) + 1)
    # 					importantFeaturesFile.write('>' + featureName + '\n')
    # 					importantFeaturesFile.write(str(feature.extract(record).seq) + '\n')
    # 					listOfImportantFeatures[featureName] = feature
    # 				else:
    # 					importantFeaturesFile.write('>' + featureName + '\n')
    # 					importantFeaturesFile.write(str(feature.extract(record).seq) + '\n')
    # 					listOfImportantFeatures[featureName] = feature

    # running blast
    if nbrRNA > 0:
        logging.info("Formatting database for blast...")
        command = "makeblastdb -in important_features.fasta -dbtype nucl"  # need to formatdb refseq first

        args = shlex.split(command)
        formatDB = Popen(args, stdout=open(os.devnull, "wb"))
        formatDB.wait()

        with open("important_features.blast.xml", "w") as blastResultFile:
            if (
                usedOwnGenBankReference == True
            ):  # using a personal genbank reference, make e-value more restrict
                command = (
                    "blastn -db important_features.fasta -query "
                    + resultFile
                    + " -outfmt 5 -num_threads 2 -word_size 8 -perc_identity "
                    + str(cutoffEquality_nucl)
                    + " -max_hsps 5 -gapextend 2 -gapopen 2 "
                    + "-dust no"
                )  # call BLAST with XML output
            else:  # using a non personal genbank reference
                command = (
                    "blastn -db important_features.fasta -query "
                    + resultFile
                    + " -outfmt 5 -num_threads 2 -word_size 8 -perc_identity "
                    + str(cutoffEquality_nucl)
                    + " -max_hsps 5 -gapextend 2 -gapopen 2 "
                    + "-dust no"
                )  # call BLAST with XML output
            args = shlex.split(command)
            blastAll = Popen(args, stdout=blastResultFile)
            blastAll.wait()
        """with open("important_features.blast2.out",'w') as blastResultFile:
			if usedOwnGenBankReference == True: #using a personal genbank reference, make e-value more restrict
				command = "blastn -db important_features.fasta -query " + resultFile + " -outfmt 6 -num_threads 2 -word_size 8 -perc_identity " + str(cutoffEquality_nucl) + " -max_hsps 5 -gapextend 2 -gapopen 2 "+ "-dust no" #call BLAST with XML output
			else: #using a non personal genbank reference
				command = "blastn -db important_features.fasta -query " + resultFile + " -outfmt 6 -num_threads 2 -word_size 8 -perc_identity " + str(cutoffEquality_nucl) + " -max_hsps 5 -gapextend 2 -gapopen 2 " + "-dust no" #call BLAST with XML output
			args = shlex.split(command)
			blastAll = Popen(args, stdout=blastResultFile)
			blastAll.wait()"""

        # SearchIO object handler and checker for best hit separation
        blastparse = SearchIO.parse(
            "important_features.blast.xml", "blast-xml"
        )  # get all queries
        """for qresult in blastparse: #in each query, let's look for a good hit
			for qhit in qresult.hits:
				for hsp in qhit.hsps: #hsp object checking, this contains the alignment info 
					featureName = qhit.id
					if alignLen >= len(targetFeature) * alignCutOff/100:
						featureFrame = hsp.hit_frame
						seqName = featureName
						alignment = Alignment(featureName, seqName, alignLen)
						alignment.refSeq = refSeq
						alignment.frame = featureFrame
						startBase = min(hsp.hit_range[0],hsp.hit_range[1])
						endBase = max(hsp.hit_range[0],hsp.hit_range[1])
						if alignLen <= len(targetFeature) * 0.98:
							queryStart = min(hsp.query_range[0],hsp.query_range[1])
							queryEnd = max(hsp.query_range[0],hsp.query_range[1])
							newEnd = endBase + (len(targetFeature) - queryEnd)
							if newEnd <= len(refSeq.seq):
								endBase = endBase + (len(targetFeature) - queryEnd)
							else:
								endBase = len(refSeq.seq)
							startBase = max(1,startBase - queryStart)
						alignment.startBase = startBase
						alignment.endBase = endBase
						alignment.seqFound = refSeq.seq[startBase:endBase]
						listOfPresentFeatures[featureName] = (listOfImportantFeatures[featureName], alignment, featureFrame == -1)
						break """

        for qresult in blastparse:  # in each query, let's look for a good hit
            for qhit in qresult.hits:
                for (
                    hsp
                ) in qhit.hsps:  # hsp object checking, this contains the alignment info
                    featureName = qhit.id
                    if float(str(hsp.ident_num) + ".00") / float(
                        str(hsp.aln_span) + ".00"
                    ) * 100 >= float(cutoffEquality_prot):
                        if featureName in listOfImportantFeatures:
                            logging.info(featureName)
                            targetFeature = listOfImportantFeatures[featureName]
                            if (
                                hsp.aln_span * 3
                                >= len(targetFeature) * alignCutOff / 100
                            ):
                                startBase = (
                                    min(hsp.query_range[0], hsp.query_range[1]) + 1
                                )
                                endBase = max(hsp.query_range[0], hsp.query_range[1])
                                alignLen = (endBase - startBase) + 1
                                if alignLen >= len(targetFeature) * alignCutOff / 100:
                                    featureFrame = hsp.hit_frame
                                    seqName = featureName
                                    alignment = Alignment(
                                        featureName, seqName, alignLen
                                    )
                                    alignment.refSeq = refSeq
                                    alignment.translationTable = organismType
                                    alignment.frame = featureFrame
                                    alignment.startBase = startBase
                                    alignment.endBase = endBase
                                    alignment.seqFound = refSeq.seq[
                                        (startBase - 1) : endBase
                                    ]
                                    listOfPresentFeatures[featureName] = (
                                        listOfImportantFeatures[qhit.id],
                                        alignment,
                                        featureFrame <= -1,
                                    )
                                    break

        shutil.copyfile("important_features.blast.xml", out_blast + "_ref.blast.xml")
        shutil.copyfile("important_features.fasta", out_blast + "_ref.fasta")
        os.remove("important_features.blast.xml")

    os.remove("important_features.fasta")

    return (
        listOfPresentFeatures,
        listOfImportantFeatures,
        listOfSplits,
        listOfCompleteGenes,
    )


def createImageOfAnnotation(sequenceObject, outputFile):
    """Creates an image of the annotation, with relative positions of features and their size"""
    try:
        from PIL import ImageFont, Image, ImageDraw
    except:
        logging.info(
            "Could not import Image or ImageDraw library, no image of result being created."
        )
        return False

    horizontalSize = 1224
    verticalSize = 250
    red = (255, 102, 102)
    green = (0, 102, 51)
    bege = (255, 178, 102)
    blue = (102, 178, 255)
    white = (255, 255, 255)
    size = (horizontalSize, verticalSize)  # size of the image to create
    im = Image.new("RGB", size, white)  # create the image
    draw = ImageDraw.Draw(im)  # create a drawing object that is
    # used to draw on the new image
    n = 1
    legenda = []

    for gbkFeature in sequenceObject.features:
        if (
            gbkFeature.type == "tRNA"
            or gbkFeature.type == "CDS"
            or gbkFeature.type == "rRNA"
            or gbkFeature.type == "D-loop"
        ):
            featureLen = (gbkFeature.location.end - gbkFeature.location.start) + 1
            featureRelativeSize = horizontalSize * featureLen / len(sequenceObject.seq)
            featureRelativeStart = (
                horizontalSize * gbkFeature.location.start / len(sequenceObject.seq)
            ) + 1

            if gbkFeature.location.strand == 1:
                if n % 2 == 0:
                    text_pos = (
                        featureRelativeStart - 1,
                        20,
                    )  # top-left position of our text
                else:
                    text_pos = (
                        featureRelativeStart - 1,
                        10,
                    )  # top-left position of our text
            else:
                if n % 2 == 0:
                    text_pos = (
                        featureRelativeStart - 1,
                        125,
                    )  # top-left position of our text
                else:
                    text_pos = (
                        featureRelativeStart - 1,
                        115,
                    )  # top-left position of our text

            for qualifier in gbkFeature.qualifiers:
                if qualifier == "product" or qualifier == "gene":
                    # get feature name
                    text = str(n)  # gbkFeature.qualifiers[qualifier]

                    if gbkFeature.qualifiers[qualifier] not in legenda:
                        legenda.append(gbkFeature.qualifiers[qualifier])

                    if gbkFeature.type == "rRNA":
                        triangleColor = red
                    elif gbkFeature.type == "tRNA":
                        triangleColor = bege
                    else:
                        triangleColor = blue

            module_dir = os.path.dirname(__file__)
            module_dir = os.path.abspath(module_dir)
            font_full_path = os.path.join(module_dir, "fonts/FreeSans.ttf")

            font = ImageFont.truetype(font_full_path, 12)

            # Now, we'll do the drawing:
            draw.text(text_pos, text, fill="black", font=font)

            if gbkFeature.location.strand == 1:
                draw.polygon(
                    [
                        (featureRelativeStart, 40),
                        (featureRelativeStart, 70),
                        (featureRelativeStart + featureRelativeSize, 55),
                    ],
                    outline=triangleColor,
                    fill=triangleColor,
                )
            else:
                draw.polygon(
                    [
                        (featureRelativeStart, 90),
                        (featureRelativeStart + featureRelativeSize, 105),
                        (featureRelativeStart + featureRelativeSize, 75),
                    ],
                    outline=triangleColor,
                    fill=triangleColor,
                )

            n += 1

    nlegenda = 0
    legendaString = ""
    linha = 1

    while draw.textsize(legendaString, font=font)[
        0
    ] < horizontalSize and nlegenda < len(legenda):
        nlegenda += 1
        if (
            draw.textsize(
                legendaString + str(nlegenda) + "-" + legenda[nlegenda - 1] + ", ",
                font=font,
            )[0]
            > horizontalSize
        ):
            draw.text(
                (0, 155 + 20 * (linha - 1)), legendaString, fill="black", font=font
            )
            linha += 1
            nlegenda -= 1
            legendaString = ""
        elif nlegenda == len(legenda):
            legendaString += str(nlegenda) + "-" + legenda[nlegenda - 1]
            draw.text(
                (0, 155 + 20 * (linha - 1)), legendaString, fill="black", font=font
            )
            break
        else:
            legendaString += str(nlegenda) + "-" + legenda[nlegenda - 1] + ", "

    draw.text(
        (horizontalSize / 2, verticalSize - 15),
        sequenceObject.name,
        fill="black",
        font=font,
    )

    del draw  # I'm done drawing so I don't need this anymore

    # now, we tell the image to save as a PNG to the
    # provided file-like object
    im.save(outputFile, "PNG")


def main():
    # Config logging when called as main.
    logging.basicConfig(
        level=0, format="%(asctime)s:%(levelname)s:%(module)s:%(message)s"
    )
    blasteVal = sys.argv[7]
    percent_equality_prot = sys.argv[8]
    percent_equality_nucl = sys.argv[9]
    genbank = sys.argv[10]
    nWalk = int(sys.argv[11])
    # This is the tool to use for tRNA search
    tRNAscan = sys.argv[12]

    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        logging.info(
            "Usage: genbank_reference fasta_file output_file organism_type(integer, default=2) alignCutOff(float, default=45) coveCutOff(7)"
        )
        logging.info("Only the first, second, and third arguments are required.")
    else:
        is_avail(["makeblastdb", "blastn", "blastx"])

        fastaReference = sys.argv[1]
        resultFile = sys.argv[2]
        out_blast = resultFile.split(".fasta")[0]
        outputFile = sys.argv[3]
        try:
            organismType = int(sys.argv[4])
            logging.info("Organism type specified: %s" % organismType)
        except:
            organismType = 2
            logging.info(
                "organism_type was not specified, assuming 2 (vertebrate mitochondria)"
            )
        try:
            alignCutOff = float(sys.argv[5])
            logging.info("alignCutOff: %s" % alignCutOff)
        except:
            alignCutOff = 45
            logging.info("alignCutOff was not specified, assuming 0.5")
        try:
            coveCutOff = int(sys.argv[6])
            logging.info("coveCutOff: %s" % coveCutOff)
        except:
            coveCutOff = 7
            logging.info("coveCutOff was not specified, assuming 7")
        x = geneCheck(
            fastaReference=fastaReference,
            resultFile=resultFile,
            cutoffEquality_prot=percent_equality_prot,
            cutoffEquality_nucl=percent_equality_nucl,
            usedOwnGenBankReference=True,
            organismType=organismType,
            alignCutOff=alignCutOff,
            out_blast=out_blast,
            blasteVal=blasteVal,
        )

        logging.info("Features found: %s" % len(x[0]))
        logging.info("Total features: %s" % len(x[1]))
        logging.info("Running tRNA annotation with " + tRNAscan)
        presentFeatures = x[0]
        assemblyCheck = tRNAscanChecker.tRNAscanCheck(
            resultFile, True, False, organismType, coveCutOff, False, False, tRNAscan
        )  # returns a Assembly object with statistics and alignment info
        tRNAs = assemblyCheck.tRNAs

        listOfFeaturesToOutput = []
        listOfFoundTRNAs = []
        for foundFeature in presentFeatures:
            thisFeatureFound = presentFeatures[foundFeature][1]
            # comparing tRNAscan-SE results with this, in case tRNAscan-SE was run
            if "trn" in thisFeatureFound.seq2.lower():
                for tRNAFound in tRNAs:
                    # down here we update the start and end positions of tRNAs found with Needle, with the
                    # results outputted by tRNAScan-SE
                    # tRNAconver = guarantees all tRNA names are in tRNA-Phe format
                    if "trna-" + tRNAFound.tRNAtype.lower() == tRNAconvert(
                        thisFeatureFound.seq2.lower()
                    ):
                        thisFeatureFound.startBase = min(
                            tRNAFound.tRNAcoordinates[0], tRNAFound.tRNAcoordinates[1]
                        )
                        thisFeatureFound.endBase = max(
                            tRNAFound.tRNAcoordinates[0], tRNAFound.tRNAcoordinates[1]
                        )
                        if tRNAFound.tRNAcoordinates[0] > tRNAFound.tRNAcoordinates[1]:
                            thisFeatureFound.frame = -1
                        else:
                            thisFeatureFound.frame = 1

                        break

                listOfFoundTRNAs.append(thisFeatureFound.seq2.lower())

            listOfFeaturesToOutput.append(thisFeatureFound)

        # if tRNAscan-SE was run, check the tRNAs it found and input them in the features to output list
        for tRNAFound in tRNAs:
            tRNAName = "trna-" + tRNAFound.tRNAtype.lower()
            if tRNAFound.tRNAintronBegin > 0:
                logging.warning(
                    "WARNING: %s was found with an intron!" % prettyRNAName(tRNAName)
                )
            if (
                tRNAName not in tRNAconvert(listOfFoundTRNAs)
                and "trna-sec" not in tRNAName
                and "trna-sup" not in tRNAName
            ):
                newTRNAStart = tRNAFound.tRNAcoordinates[0]
                newTRNAEnd = tRNAFound.tRNAcoordinates[1]
                newTRNALen = max(newTRNAStart, newTRNAEnd) - min(
                    newTRNAStart, newTRNAEnd
                )
                newTRNA = Alignment(tRNAName, prettyRNAName(tRNAName), newTRNALen)
                newTRNA.startBase = min(newTRNAStart, newTRNAEnd) + 1
                newTRNA.endBase = max(newTRNAStart, newTRNAEnd)
                thisFeatureFound = newTRNA

                if newTRNAStart > newTRNAEnd:
                    newTRNA.frame = -1
                else:
                    newTRNA.frame = 1

                presentFeatures[prettyRNAName(tRNAName)] = (
                    False,
                    thisFeatureFound,
                    False,
                )

                listOfFeaturesToOutput.append(thisFeatureFound)

        listOfFeaturesToOutput.sort()
        logging.info(
            "Total features found after " + str(tRNAscan) + ": ",
            len(listOfFeaturesToOutput),
        )

        finalResults = genbankOutput.genbankOutput(
            outputFile, resultFile, listOfFeaturesToOutput, False, 900, nWalk
        )

        with open(outputFile, "w") as outputResult:
            count = SeqIO.write(finalResults, outputResult, "genbank")
            createImageOfAnnotation(finalResults, outputFile.split(".gb")[0] + ".png")

        outputFile = outputFile.split(".gb")[0] + "_raw.gff"
        outputFile = open(outputFile, "w")
        seq = SeqIO.read(open(resultFile, "r"), "fasta")
        seq_name = seq.name

        genes = {}
        for gbkFeature in finalResults.features:
            for qualifier in gbkFeature.qualifiers:
                if qualifier == "product" or qualifier == "gene":
                    if gbkFeature.location.strand == 1:
                        direction = "+"
                    if gbkFeature.location.strand == -1:
                        direction = "-"
                    if not gbkFeature.qualifiers[qualifier] in genes:
                        outputFile.write(
                            seq_name
                            + "\t"
                            + "mitofinder"
                            + "\t"
                            + str(gbkFeature.type)
                            + "\t"
                            + str(gbkFeature.location.start + 1)
                            + "\t"
                            + str(gbkFeature.location.end)
                            + "\t"
                            + "."
                            + "\t"
                            + direction
                            + "\t"
                            + "0"
                            + "\t"
                            + str(gbkFeature.qualifiers[qualifier])
                            + "\n"
                        )
                        genes[gbkFeature.qualifiers[qualifier]] = gbkFeature.qualifiers[
                            qualifier
                        ]
        outputFile.close()
        """
		if ('TRNF' in presentFeatures) or ('tRNA-Phe' in presentFeatures) or ('trnf' in presentFeatures) or ('trnF' in presentFeatures):
			logging.info('Creating ordered genbank file (with tRNA-Phe at the start)...')
			resultOrderedGbFile = outputFile.replace('.gb','') + '.ordered.gb'

			if 'TRNF' in presentFeatures:
				lookForPhe = 'TRNF'
			elif 'tRNA-Phe' in presentFeatures:
				lookForPhe = 'tRNA-Phe'
			elif 'trnf' in presentFeatures:
				lookForPhe = 'trnf'
			elif 'trnF' in presentFeatures:
				lookForPhe = 'trnF'

			pheAlignment = presentFeatures[lookForPhe][1]
			pheStart = pheAlignment.startBase
			orderedFinalResults = finalResults[pheStart:] + finalResults[0:pheStart]

			with open(resultOrderedGbFile, "w") as outputResult: #create the file!
				count = SeqIO.write(orderedFinalResults, outputResult, "genbank")
				count = SeqIO.write(orderedFinalResults, resultOrderedGbFile.replace('.gb','.fasta'), "fasta")
				createImageOfAnnotation(orderedFinalResults, 'orderedResult.png')
			
		else:
			logging.info("Since tRNA-Phe couldn't be found, ordered genbank file wasn't created.")
		"""


if __name__ == "__main__":
    main()
