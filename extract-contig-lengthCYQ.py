#this script was used to calculate the length of contigs that list in certain files.
from Bio import SeqIO
inputfile=input("Enter the putative virus contigs (>5 viral proteins) file: ")
assemblyfile=input("Enter the corresponding assembed files: ")
outputfile=open("length_of_"+inputfile,'w')
outputfile2=open("over_5000bp_"+inputfile,'w')
a=[]
for line in open(inputfile,'r'):
    a.append(str(line).strip())

for record in SeqIO.parse(open(assemblyfile,'r'),'fasta'):
    if str(record.id).strip() in a:
        outputfile.write(str(record.id).strip()+'\t'+str(len(str(record.seq).strip()))+'\n')
        if float(len(str(record.seq))) >= 5000:
            outputfile2.write(str(record.id).strip()+'\t'+str(len(str(record.seq).strip()))+'\n')
        else:
            continue
    else:
        continue

outputfile.close()
outputfile2.close()
input("Press <enter> to close the window!")
