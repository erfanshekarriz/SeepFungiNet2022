"""
Scans multiple assembled transcriptome files in fasta or fa format within a single input directory
and finds highly probable hydrophobins genes.
"""

#! /usr/local/bin/python3.10

import os
import sys
import re


inputdir = input("Enter the input file directory with all fasta files: ")
seqtype = input("Enter whether your sequence type (RNA or AA): ")
seqtype = seqtype.lower()

# Read input file fasta file function


def fasta2dict(fastafile):
    """Reads a fasta file and returns a dictionary"""
    try:
        with open(fastafile, "r", errors="ignore", encoding="UTF-8") as inputfile:
            content = inputfile.read()
            contentsplit = content.split(">")

            cleanlist = []
            for seql in contentsplit:
                seqr = repr(seql)
                clean = seqr.replace("\\n", "*NEWLINE*", 1)
                clean2 = clean.replace("\\n", "")
                clean3 = clean2.replace("'", "")
                cleanlist.append(clean3)

            seq_dict = {}
            for seqc in cleanlist:
                seqto = seqc
                seqsplit = seqto.split("*NEWLINE*")
                if len(seqsplit) > 1:
                    seqkey = ">" + seqsplit[0]
                    seq_dict[seqkey] = seqsplit[1]

            return seq_dict

    except FileNotFoundError:
        print("Could not open/read file: ", fastafile)


# check to see if file has already been procssessed
# RNA CHECK
if seqtype == "rna":
    # make new output directory
    os.chdir(inputdir)
    outputdir = os.path.join(inputdir, "motifchecked_RNA")
    os.mkdir(outputdir)

    for filename in os.listdir(inputdir):
        os.chdir(inputdir)
        print(filename)
        root, ext = os.path.splitext(filename)
        filecheck = outputdir + "/motifchecked_" + root + ext

        # check to see if file is fasta or fa format and if file already exists
        # if requirements are met then continue with processing
        if filename == "motifchecked":
            continue

        elif ext in (".fa", ".fasta", ".fna", ".fas") and not os.path.isfile(filecheck):

            # read fasta file
            finalseqDict = {}
            SeqDict = fasta2dict(filename)

            # use regular expressions to match hydrophobin cystein patterns
            # --C--CC--C--C--CC--C--
            STRPAT = (
                "^.+?(((TGC|TGT)).+?((TGC|TGT)(TGC|TGT)).+?((TGC|TGT)).+?"
                "((TGC|TGT)).+?((TGC|TGT)(TGC|TGT)).+?((TGC|TGT))).+$"
            )
            cyspattern = re.compile(STRPAT, re.IGNORECASE)

            for seqID, seq in SeqDict.items():
                cysMatch = re.search(cyspattern, seq)

                # make sure that the cysteine patterns are on the same reading frame
                if cysMatch is not None:
                    # The first cysteine is the reference frame
                    ref_frame = cysMatch.start(2)

                    matchgrpindx = [
                        cysMatch.start(4),
                        cysMatch.start(7),
                        cysMatch.start(9),
                        cysMatch.start(11),
                        cysMatch.start(14),
                    ]

                    readfraemlst = [((indx - ref_frame) % 3) for indx in matchgrpindx]

                    if len(set(readfraemlst)) == 1 and readfraemlst[0] == 0:
                        finalseqDict[seqID] = seq

                else:
                    continue

            os.chdir(outputdir)
            with open(
                ("motifchecked_" + root + ext), "w", encoding="utf-8"
            ) as outputfile:
                for seqID, seq in finalseqDict.items():
                    outputfile.write(seqID + "\n" + seq + "\n")
        else:
            print()
            print(
                "\n"
                + os.path.basename(filename)
                + " has already been scanned or is not in the correct .fa or .fasta format..."
            )
            print("Continuing to next file!\n")
            continue

# AMINO ACID CHECK
elif seqtype == "aa":
    # make new output directory
    os.chdir(inputdir)
    outputdir = os.path.join(inputdir, "motifchecked_Prot")
    os.mkdir(outputdir)

    for filename in os.listdir(inputdir):
        os.chdir(inputdir)
        print(filename)
        root, ext = os.path.splitext(filename)
        filecheck = outputdir + "/motifchecked_" + root + ext

        # check to see if file is fasta or fa format and if file already exists
        # if requirements are met then continue with processing
        if filename == "motifchecked":
            continue

        elif ext in (".fa", ".fasta", ".fna", ".fas") and not os.path.isfile(filecheck):

            # read fasta file
            finalseqDict = {}
            SeqDict = fasta2dict(filename)

            # use regular expressions to match hydrophobin cystein patterns
            # --C--CC--C--C--CC--C--
            STRPAT = "^.+?C.+?CC.+?C.+?C.+?CC.+?C.+$"
            cyspattern = re.compile(STRPAT, re.IGNORECASE)
            for seqID, seq in SeqDict.items():
                cysMatch = re.search(cyspattern, seq)
                cyscount = seq.count("C")
                # print(cyscount)


                # make sure that the cysteine patterns are on the same reading frame
                if (cysMatch is not None): #and (8 <= cyscount <= 10):
                    # The first cysteine is the reference frame
                    # print(seq)
                    # print(cyscount)
                    finalseqDict[seqID] = seq

                else:
                    continue

            os.chdir(outputdir)
            with open(
                ("motifchecked_" + root + ext), "w", encoding="utf-8"
            ) as outputfile:
                for seqID, seq in finalseqDict.items():
                    outputfile.write(seqID + "\n" + seq + "\n")
        else:
            print()
            print(
                "\n"
                + os.path.basename(filename)
                + " has already been scanned or is not in the correct .fa or .fasta format..."
            )
            print("Continuing to next file!\n")
            continue

else:
    sys.exit('Invalid sequence type. Please either choose "RNA" or "AminoAcid"')
