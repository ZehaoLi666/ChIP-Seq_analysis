fa="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
gff="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"
seqs={};name=None;buf=[]
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name]="".join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
seqs[name]="".join(buf)
S=seqs["Pf3D7_10_v3"].upper(); L=len(S)
print("N count chr10:", S.count("N"))
# N runs
runs=[];i=0
while i<L:
    if S[i]=="N":
        j=i
        while j<L and S[j]=="N": j+=1
        if j-i>=50: runs.append((i,j,j-i))
        i=j
    else: i+=1
print("N runs >=50bp:", runs[:10])

# genes on chr10, find largest intergenic gaps
g=[]
for line in open(gff):
    if line.startswith("#"): continue
    f=line.split("\t")
    if len(f)>4 and f[2]=="gene" and f[0]=="Pf3D7_10_v3":
        g.append((int(f[3])-1,int(f[4])))
g.sort()
gaps=[]
prev=0
for s,e in g:
    if s>prev: gaps.append((prev,s,s-prev))
    prev=max(prev,e)
gaps.sort(key=lambda x:-x[2])
print("\ntop 8 intergenic gaps on chr10 (start,end,size):")
for gg in gaps[:8]:
    mid=(gg[0]+gg[1])//2
    w=S[max(0,mid-1000):mid+1000]
    at=(w.count("A")+w.count("T"))/len(w)
    print(f"  {gg[0]}-{gg[1]} size={gg[2]} mid={mid} AT_2kb_at_mid={at:.4f}")

print("\ngenes overlapping 925-950 kb:")
for s,e in g:
    if e>925000 and s<950000: print(f"  {s}-{e}")
