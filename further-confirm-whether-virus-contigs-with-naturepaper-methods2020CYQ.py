#This manuscript was used to confirm the taxon of contigs that have more than 5 viral ORFs
from Bio import SeqIO

inputfile1=input("Enter the filename in which contains contigs ids (>5kb&>=5viral ORFs): ")
inputfile2=input("Enter the filename in which contains the putative contigs' Pfam hmmsearch results: ")
inputfile3=input("Enter the filename in which contains the putative contigs' KEGG diamond blast results: ")
inputfile4=input("Enter the origanl assembled sequence at aa format: ")#the original assemly file with aa format
inputfile5=input("Enter the filename in which contains all ORFs' hmmsearch results,which search against the nature paper database: ")

#------------------------------************************--------------------------------------******************************
a=[]
#count the number of ORFs of each assembled contig.
for record in SeqIO.parse(open(inputfile4,'r'),'fasta'):
    content=str(record.id).strip().split('_')
    ctid1=str(content[0])#contig ids of all assemlbed contigs
    a.append(ctid1)
b=list(set(a))
b.sort() #store all putative viral contigs

c={}#store the orf number of each contig
for item in b:
    c[str(item)]=str(a.count(item)) #c is the library that contains information of contig id-orf number
temp=open('contig-orf-number.txt'+inputfile4,'w')
for key,value in c.items():
    temp.write(str(key)+'\t'+str(value)+'\n')
temp.close()

#------------------------------************************--------------------------------------******************************
#collect all priliminary contig ids that have 5 viral hits and > 5 kbp length.

CT=[] #store all putative viral contigs based on the preliminary analysis
for line in open(inputfile1,'r'):
    ctid = str(line).strip().split('\t')[0]
    CT.append(ctid)

#------------------------------************************--------------------------------------******************************
#ORF with viral hits, Pfam hits, and KEGG hits; find the overlapped ORFs.
viral=[]
pfam=[]
KEGG=[]
for line in open(inputfile5,'r'):
    if line.startswith('#'):
        continue
    else:
        orfid1=str(line).strip().split(' ')[0]
        viral.append(orfid1)
tpviral=list(set(viral)) #list includes all ORFs with viral hits

for line in open(inputfile2,'r'):
    if line.startswith('#'):
        continue
    else:
        orfid2=str(line).strip().split(' ')[0]
        pfam.append(orfid2)
tppfam=list(set(pfam)) #list includes all ORFs with Pfam hits.

for line in open(inputfile3,'r'):
    orfid3=str(line).strip().split('\t')[0]
    KEGG.append(orfid3)
tpKEGG=list(set(KEGG)) #list includes all ORFs with KEGG hits

olpfam=list((set(tppfam).union(set(tpviral)))^(set(tppfam)^set(tpviral))) #find overlapped orf ids (pfam)
olKEGG=list((set(tpKEGG).union(set(tpviral)))^(set(tpKEGG)^set(tpviral))) #find overlapped orf ids (KEGG)

#------------------------------************************--------------------------------------******************************
#find the real viral ORFs after removing the overlapped orfs, and store the real viral ORFs number of putative virus contigs
overlap=[]
for item in olpfam:
    overlap.append(str(item))
for item in olKEGG:
    overlap.append(str(item))
mightrm=list(set(overlap)) #ORFs might be removed

for item in mightrm:
    tpviral.remove(str(item))

virallib={} #store the real viral orf number of purified virus contigs

viralct=[] #contigs contian real viral orfs (reduntant)
for item in tpviral:
    content=str(item).strip().split('_')
    ctid2=str(content[0])
    viralct.append(str(ctid2))
virct=list(set(viralct)) #contigs contian real viral orfs (nonreduntant)

for item in CT:
    if str(item) in virct:
        virallib[str(item)]=str(viralct.count(str(item)))
    else: #contigs may do not have real viral ORFs
        virallib[str(item)]='0'

#------------------------------************************--------------------------------------******************************
#calculate the Pfam hits of each putative viral contig

Pfamlib={}#store the number of Pfam terms of each putative contig
Pfamorf=[] #store all ORF ids with Pfam hits (redundanct)
Pfamct=[] #store all contig ids carring ORFs with Pfam hits (redundant)

for line in open(inputfile2,'r'):
    if line.startswith('#'):
        continue
    else:
        orfid=str(line).strip().split(' ')[0]
        content=orfid.split('_')
        ctid3=str(content[0])
        Pfamorf.append(str(orfid))
        Pfamct.append(str(ctid3))
Pfamorfuq=list(set(Pfamorf)) #store all ORF ids with Pfam hits (non-redundanct)
Pfamctuq=list(set(Pfamct)) #store all contig ids carring ORFs with Pfam hits (non-redundant)

tmp1=[]
for item in Pfamorfuq:
    content=str(item).strip().split('_')
    ctid4=str(content[0])
    tmp1.append(str(ctid4))
Pfamctuq.sort()
for item in CT:
    if item in Pfamctuq:
        Pfamlib[str(item)]=tmp1.count(str(item))
    else:#contigs may do not have orf matched with pfam database
        Pfamlib[str(item)]='0'

#------------------------------************************--------------------------------------******************************
#calculate the KEGG hits of each putative viral contigs

KEGGlib={} #store the number of KEGG terms of each putative contig
KEGGorf=[] #store all ORF ids with KEGG hits (redundanct)
KEGGct=[] #store all contig ids with KEGG hits (redundanct)

for line in open(inputfile3,'r'):
    orfid=str(line).strip().split('\t')[0]
    content=str(orfid).strip().split('_')
    ctid5=str(content[0])
    KEGGct.append(str(ctid5))
    KEGGorf.append(str(orfid))
KEGGorfuq=list(set(KEGGorf)) #store all ORF ids with KEGG hits (non-redundanct)
KEGGctuq=list(set(KEGGct)) #store all contig ids with KEGG hits (non-redundanct)

tmp2=[]
for item in KEGGorfuq:
    content=str(item).strip().split('_')
    ctid6=str(content[0])
    tmp2.append(str(ctid6))
KEGGctuq.sort()
for item in CT:
    if item in KEGGctuq:
        KEGGlib[str(item)]=tmp2.count(str(item))
    else:#contigs may do not have orf matched with pfam database
        KEGGlib[str(item)]='0'

#------------------------------************************--------------------------------------******************************
# 1st criterion
outputfile1=open('1st-step_Meet_with_Pfam_cutoff'+inputfile1,'w')
outputfile1.write("contig-id"+"\t"+"viral-number"+"\t"+"KO-proportation"+"\t"+"Pfam-proportation"+'\n')
for item in CT:
    if float(virallib[str(item)])/float(c[str(item)])>=0.1:
        if float(KEGGlib[str(item)])/float(c[str(item)])<=0.2:
            if float(Pfamlib[str(item)])/float(c[str(item)])<=0.4:
                print(float(KEGGlib[str(item)])/float(c[str(item)]))
                print(float(Pfamlib[str(item)])/float(c[str(item)]))
                outputfile1.write(str(item)+'\t'+str(virallib[str(item)])+'\t'+str(float(KEGGlib[str(item)])/float(c[str(item)]))+'\t'+str(float(Pfamlib[str(item)])/float(c[str(item)]))+'\n')
            else:
                continue
        else:
            continue
    else:
        continue
outputfile1.close()

#------------------------------************************--------------------------------------******************************
# 2nd and 3rd criteria
outputfile2=open('2ndstep_Meet_with_cutoff'+inputfile1,'w')#check whether the viral ORFs>=Pfam
outputfile3=open('3rdstep_Meet_with_cutoff'+inputfile1,'w')#check whether the vrial ORFs>60%ORFs
outputfile2.write("contig-id"+"\t"+"viral-ORFs"+"\t"+"Pfam-ORFs"+"\n")
outputfile3.write("contig-id"+"\t"+"viral-ORFs/total-ORFs"+"\n")
for item in CT:
    if float(virallib[str(item)]) >= float(Pfamlib[str(item)]):
        outputfile2.write(str(item)+'\t'+str(float(virallib[str(item)]))+'\t'+str(float(Pfamlib[str(item)]))+'\n')
    if float(virallib[str(item)])/float(c[str(item)])>=0.6:
        outputfile3.write(str(item)+'\t'+str(float(virallib[str(item)])/float(c[str(item)]))+'\n')
    else:
        continue
outputfile2.close()
outputfile3.close()
input("Press <enter> to close this window.")
