fa="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
seqs={};name=None;buf=[]
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name]="".join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
seqs[name]="".join(buf)
S=seqs["Pf3D7_10_v3"].upper(); L=len(S); W=2000
cands=[]
for i in range(60000, L-W-60000, 200):
    w=S[i:i+W]; at=(w.count("A")+w.count("T"))/W
    cands.append((at,i))
cands.sort(reverse=True)
# report top non-overlapping peaks
kept=[]
for at,i in cands:
    if all(abs(i-j)>10000 for _,j in kept):
        kept.append((at,i))
    if len(kept)>=6: break
print("chr10 top AT windows (2 kb):")
for at,i in kept: print(f"  {i}-{i+W}\tAT={at:.4f}")
