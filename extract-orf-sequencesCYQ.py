from Bio import SeqIO

inputfile = input("Enter the orf faa file: ")
idfile = input("Enter the input file which contains putative viral contigs and contig length: ")
outputfile = open("Extracted-orfs-"+str(inputfile),"w")

a =[] #store all putative viral contigs
for line in open(idfile,"r"):
    a.append(str(line).strip().split("\t")[0])
i=0
for record in SeqIO.parse(open(inputfile,"r"),"fasta"):
    orfid = str(record.id)
    ctid = str(orfid).strip().split("_")[0]
    if str(ctid) in a:
        i+=1
        outputfile.write(">"+str(orfid)+'\n')
        outputfile.write(str(record.seq).strip()+'\n')
    else:
        continue

outputfile.close()
print("There are ",i," orfs extraced from",str(inputfile))

input("Press <enter> to close this window.")
