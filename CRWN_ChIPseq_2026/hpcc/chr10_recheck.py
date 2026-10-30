fa="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
seqs={};name=None;buf=[]
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name]="".join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
seqs[name]="".join(buf)
S=seqs["Pf3D7_10_v3"].upper(); L=len(S)
print("chr10 length", L)
print("=== 2-kb windows, 900-980 kb ===")
for i in range(900000, 980000, 2000):
    w=S[i:i+2000]; at=(w.count("A")+w.count("T"))/2000
    print(f"  {i}  AT={at:.4f}")
print("=== best 2-kb window genome-position scan on chr10, excluding 60kb ends ===")
cands=[]
for i in range(60000, L-2000-60000, 500):
    w=S[i:i+2000]; cands.append(((w.count("A")+w.count("T"))/2000, i))
cands.sort(reverse=True)
kept=[]
for at,i in cands:
    if all(abs(i-j)>15000 for _,j in kept): kept.append((at,i))
    if len(kept)>=8: break
for at,i in kept: print(f"  {i}-{i+2000}  AT={at:.4f}")
