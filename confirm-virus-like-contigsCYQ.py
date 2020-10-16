inputfile=input("Enter the hmmsearch results that only contain the orf ids: ")
outputfile1=open("putative_virus_contigs_"+inputfile,'w')
outputfile2=open("contigs_with_less_than_5_viral_hits_"+inputfile,'w')
a=[]
b=[]
for line in open(inputfile,'r'):
    cgid=str(line).strip().split('_')[0]
    a.append(cgid)
    orfid=str(line).strip()
    b.append(orfid)
uqorf = list (set(b))

uncg = list (set(a))

c=[]
for item in uqorf:
    cgid2=str(item).strip().split('_')[0]
    c.append(str(cgid2))
for i in uncg:
    if float(c.count(str(i)))>= 5:
        outputfile1.write(str(i).strip()+'\n')
    else:
        outputfile2.write(str(i).strip()+'\n')
outputfile1.close()
outputfile2.close()
input("Press <enter> to close this window.")
