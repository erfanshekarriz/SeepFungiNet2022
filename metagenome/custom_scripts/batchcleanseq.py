import os

inputdir = input("Input the fasta file directory: ")
filetype = ".fa"

SeqDictC = {}

os.chdir(inputdir)
for filename in os.listdir(inputdir):
    root, ext = os.path.splitext(filename)
    if ext in (".fa", ".fasta", ".fna", ".fas"):
        with open(filename,"r") as fi:
            content = fi.read()
            contentclean = content.replace("-", "")
            contentsplit = contentclean.split(">")

        cleanlist = []
        for seq in contentsplit:
            seqr = repr(seq)
            clean = seqr.replace("\\n","*NEWLINE*",1)
            clean2 = clean.replace("\\n", "")
            clean3 = clean2.replace("'","")
            cleanlist.append(clean3)

        for seq in cleanlist:
            seqto = seq
            seqsplit = seqto.split("*NEWLINE*")
            if len(seqsplit)>1:
                seqID = ">" + seqsplit[0]
                SeqDictC[seqID]= seqsplit[1]
        #print(SeqDictC)

        with open(("clean_" + filename), "w+") as f:
            for seqID, seqz in SeqDictC.items():
                seqIDnew = (seqID.replace("_", " ")).replace(".", " ").replace("|", " ")
                seqIDnew = (seqIDnew.replace("[", " ")).replace(":", " ").replace("+", " ")
                seqIDnew = (seqIDnew.replace("]", " ")).replace("-", " ").replace("_", " ")
                seqIDnew = (seqIDnew.replace("(", " ")).replace(")", " ")

                f.write(seqIDnew + "\n" + seqz + "\n")
        print ('Done with ' + filename + "!\n")
