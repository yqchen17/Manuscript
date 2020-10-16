inputfile1=input("Enter the 1ststep filtered file name: ")
inputfile2=input("Enter the 2ndstep filtered file name: ")
inputfile3=input("Enter the 3rdstep filtered file name: ")
inputfile4=input("Enter the contigs identified by EM: ")
inputfile5=input("Enter the WWTP ID：")
outputfile=open("final_viral_contig_only_nature_based"+inputfile5+".txt",'w')

first=[]
second=[]
third=[]
delete=[]

for line in open(inputfile1,'r'):
    if "contig-id" in line:
        continue
    else:
        first.append(str(line).strip().split('\t')[0])
for line in open(inputfile2,'r'):
    if "contig-id" in line:
        continue
    else:
        second.append(str(line).strip().split('\t')[0])
for line in open(inputfile3,'r'):
    if "contig-id" in line:
        continue
    else:
        third.append(str(line).strip().split('\t')[0])
for line in open(inputfile4,'r'):
    delete.append(str(line).strip().split('\t')[0])
t1=list(set(first).union(set(second)))
t2=list(set(t1).union(set(third)))
final=list(set(t2).union(set(delete)))
t2.sort(key=int)
inter=list(set(t2).intersection(set(delete)))
print("id count of nature methods: "+str(len(t2)))
print("id count of EM method: "+str(len(delete)))
print("id count of combined: "+str(len(final)))
print("id count of shared: "+str(len(inter)))

for item in t2:
    outputfile.write(str(item).strip()+'\n')
outputfile.close()
input("Press <enter> to close this window.")