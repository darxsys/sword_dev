#!/usr/bin/python

import os, sys, getopt

def main():

    options = "i:j:f:n:h"
    longOptions = ["help", "input=", "database=", "false-positive-rate=", "name="]

    inputPath = None
    databasePath = None
    fpRate = 10000
    name = "program"

    try:
        opts, args = getopt.getopt(sys.argv[1:], options, longOptions)
    except getopt.GetoptError as err:
        print(str(err))
        help()
        sys.exit()

    for o, a in opts:
        if o in ("-i", "--input"):
            inputPath = a
        elif o in ("-j", "--database"):
            databasePath = a
        elif o in ("-f", "false-positive-rate"):
            fpRate = int(a)
        elif o in ("-n", "name"):
            name = a
        elif o in ("-h", "-help"):
            help()
            sys.exit()

    if inputPath is None:
        error("missing option: -i <file>")

    if not os.path.isfile(inputPath):
        error("non-existent file: -i {}".format(inputPath))

    if databasePath is None:
        error("missing option: -j <file>")

    if not os.path.isfile(databasePath):
        error("non-existent file: -j {}".format(databasePath))

    if fpRate < 0:
        error("invalid fp rate")

    printRocData(processInput(inputPath, processDatabase(databasePath)), fpRate, name)
    
#******************************************************************************
#******************************************************************************

def processDatabase(databasePath):
    database = open(databasePath).read().split('>')

    db = {}

    for seq in database:
        if not seq:
            continue

        lines = seq.split("\n")
        seqName = lines[0]
        superFamily = ".".join(seqName.split(" ", 2)[1].split(".")[0:-1])

        db[seqName.split(" ", 1)[0]] = superFamily

    return db

def processInput(inputPath, superFamilies):
    input_ = open(inputPath).read().split("\n")

    alignments = []

    for i in range(len(input_)):
        if not input_[i] or input_[i][0] == '#':
            continue

        if "[WARNING:" in input_[i] or "Using:" in input_[i] or "Query id" in input_[i]:
            continue;

        if "CONVERGED!" in input_[i]:
            continue;

        queryName = input_[i].split(None, 1)[0]
        targetName = input_[i].split(None, 2)[1]

        evalue = float(input_[i].split()[-2:-1][0])

        alignments.append([evalue, superFamilies[queryName] == superFamilies[targetName]])

    return alignments

def printRocData(alignments, maxFpRate, name):
    alignments.sort(key=lambda x: x[0])

    fpRate = 0
    tpRate = 0

    for align in alignments:
        if fpRate == maxFpRate:
            break

        if align[1] == True:
            tpRate += 1
        else:
            fpRate += 1

        print name + ", " + str(tpRate) + ", " + str(fpRate)


#******************************************************************************
#******************************************************************************

def error(message):
    print("[ERROR] {}".format(message))
    sys.exit()

def help():
    print(
    "usage: python rocData.py [arguments ...]\n"
    "arguments:\n"
    "    -i, --input <file>\n"
    "        (required)\n"
    "        input file in bm8 format (blast m8 tabular output format)\n"
    "    -j, --database <file>\n"
    "        (required)\n"
    "        input database file in fasta format\n"
    "    -f, --false-positive-rate <int>\n"
    "        default: 10000\n"
    "        false positive rate\n"
    "    -n, --name <string>\n"
    "        default: program\n"
    "        name of executable that produced the input file\n"
    "    -h, --help\n"
    "        prints out the help")

#******************************************************************************
#******************************************************************************
if __name__ == "__main__":
    main()
