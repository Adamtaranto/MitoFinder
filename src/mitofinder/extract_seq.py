import sys
import os.path
import shlex
import subprocess


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


allCDS = open(sys.argv[1])

dico = {}
for name, seq in read_fasta(allCDS):
    gene = name.split("@")[1]
    sp = name.split("@")[0]
    if dico.has_key(gene):
        if len(seq) > len(dico.get(gene)):
            dico[gene] = seq
    else:
        dico[gene] = seq

for cle, valeur in dico.items():
    fout = open(cle + "_all_sp.fasta", "a")
    fout.write(sp + "\n" + valeur + "\n")
