#!/usr/bin/python

import os, sys, getopt

def main():

    options = "i:j:h"
    longOptions = ["help", "first=", "second="]

    firstPath = None
    secondPath = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], options, longOptions)
    except getopt.GetoptError as err:
        print(str(err))
        help()
        sys.exit()

    for o, a in opts:
        if o in ("-i", "--first"):
            firstPath = a
        elif o in ("-j", "--second"):
            secondPath = a
        elif o in ("-h", "-help"):
            help()
            sys.exit()

    if firstPath is None:
        error("missing option: -i <file>")

    if not os.path.isfile(firstPath):
        error("non-existent file: -i {}".format(firstPath))

    if secondPath is None:
        error("missing option: -j <file>")

    if not os.path.isfile(secondPath):
        error("non-existent file: -j {}".format(secondPath))

    alignments1 = processInput(firstPath)
    alignments2 = processInput(secondPath)

    print firstPath + " vs " + secondPath

    queryLen = len(alignments1) if len(alignments1) > len(alignments2) else len(alignments2)
    print "QueryLen = " + str(queryLen)

    printEvalueData(alignments1, alignments2)
    printRankData(alignments1, alignments2)
    printAlignLenData(alignments1, alignments2)
    
#******************************************************************************
#******************************************************************************

def processInput(inputPath):
    input_ = open(inputPath).read().split("\n")

    alignments = {}
    rank = 1

    sameAligments = []

    for i in range(len(input_)):
        if not input_[i] or input_[i][0] == '#':
            continue

        if "[WARNING:" in input_[i] or "Using:" in input_[i] or "Query id" in input_[i]:
            continue;

        if "CONVERGED!" in input_[i]:
            continue;

        queryName = input_[i].split('\t', 1)[0]

        if not queryName in alignments:
            alignments[queryName] = []
            sameAligments = []
            rank = 1

        targetName = input_[i].split('\t', 2)[1]

        if targetName in sameAligments:
            continue

        sameAligments.append(targetName)

        alignLen = int(input_[i].split("\t", 4)[3])
        evalue = float(input_[i].split('\t')[-2:-1][0])

        alignments[queryName].append([targetName, alignLen, evalue, rank])

        rank += 1

    return alignments

def printEvalueData(alignments1, alignments2):
    evalues = [1.0e-250, 1.0e-100, 1.0e-75, 1.0e-50, 1.0e-25, 1.0e-10, 1.0e-7, 1.0e-4, 1.0e-3, 1.0e-2, 1, 10]

    print "Evalue, Equal, %"

    for evalue in evalues:
        equal = 0
        totalLen = 0

        for queryName in alignments1:
            temp1 = []
            temp2 = []

            for align in alignments1[queryName]:
                if align[2] > evalue:
                    break
                temp1.append(align[0])

            if queryName in alignments2:
                for align in alignments2[queryName]:
                    if align[2] > evalue:
                        break
                    temp2.append(align[0])

            for t in temp1:
                if t in temp2:
                    equal += 1

            totalLen += len(temp1)

        total = 0 if totalLen == 0 else float(equal) / float(totalLen)
        print str(evalue) + ", " + str(equal) + ", " + str(total)

    print ""

def printRankData(alignments1, alignments2):
    ranks = [10, 25, 50, 75, 100, 150, 200, 300, 400, 500]

    print "Rank, Equal, %"

    for rank in ranks:
        equal = 0
        totalLen = 0

        for queryName in alignments1:
            temp1 = []
            temp2 = []

            for align in alignments1[queryName]:
                if align[3] > rank:
                    break
                temp1.append(align[0])

            if queryName in alignments2:
                for align in alignments2[queryName]:
                    if align[3] > rank:
                        break
                    temp2.append(align[0])

            for t in temp1:
                if t in temp2:
                    equal += 1

            totalLen += len(temp1)

        total = 0 if totalLen == 0 else float(equal) / float(totalLen)
        print str(rank) + ", " + str(equal) + ", " + str(total)

    print ""

def printAlignLenData(alignments1, alignments2):
    print "AlignLen Equal, %"

    equal = 0
    totalLen = 0

    for queryName in alignments1:
        temp1 = {}
        temp2 = {}

        for align in alignments1[queryName]:
            temp1[align[0]] = align[1]

        if queryName in alignments2:
            for align in alignments2[queryName]:
                temp2[align[0]] = align[1]            

        for key in temp2:
            if key in temp1:
                if abs(temp1[key] - temp2[key]) < 5:
                    equal += 1
                totalLen += 1

    total = 0 if totalLen == 0 else float(equal) / float(totalLen)
    print str(equal) + ", " + str(total)

#******************************************************************************
#******************************************************************************

def error(message):
    print("[ERROR] {}".format(message))
    sys.exit()

def help():
    print(
    "usage: python simData.py [arguments ...]\n"
    "arguments:\n"
    "    -i, --first <file>\n"
    "        (required)\n"
    "        first input file in bm8 format (blast m8 tabular output format)\n"
    "    -j, --second <file>\n"
    "        (required)\n"
    "        second input file in bm8 format (blast m8 tabular output format)\n"
    "    -h, --help\n"
    "        prints out the help")

#******************************************************************************
#******************************************************************************

if __name__ == "__main__":
    main()
