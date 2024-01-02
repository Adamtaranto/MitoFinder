import sys
import os.path

f = open(sys.argv[1])
g = open(sys.argv[1])
evalmin = sys.argv[2]

# check species are present only once
dico = {}
dicof = {}
for line in f:
    trans = line.split("	")[1]
    trans = trans.split("	")[0]
    valeur = trans
    dico[trans] = 0
    dicof[valeur] = str("0 ()")


for cle, valeur in dico.items():
    for line in open(sys.argv[1]):
        if cle in line:
            blast = line.split("	")[2]
            blast = blast.split("	")[0]
            cibl = line.split("	")[0]
            blast = float(blast)
            eVal = line.split("	")[10]
            eVal = eVal.split("	")[0]
            test = dicof.get(cle)
            test = test.split(" ")[0]
            if float(blast) > float(test) and float(eVal) > float(evalmin):
                dicof[cle] = str(float(blast)) + " (" + cibl + ")"
for cle, valeur in dicof.items():
    print(cle + " === " + str(valeur))
