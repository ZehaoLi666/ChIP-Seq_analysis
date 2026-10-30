fa="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
seqs={};name=None;buf=[]
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name]="".join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
seqs[name]="".join(buf)
# 4-kb window (matches published PfCENH3 domain size), all chroms, for comparison
print("chrom\tbest4kb_start\tbest4kb_AT\trank2_start\trank2_AT")
for ch in sorted(seqs):
    if not ch.startswith("Pf3D7_") or "API" in ch: continue
    S=seqs[ch].upper(); L=len(S)
    if L<200000: continue
    cands=[]
    for i in range(60000, L-4000-60000, 500):
        w=S[i:i+4000]; cands.append(((w.count("A")+w.count("T"))/4000, i))
    cands.sort(reverse=True)
    kept=[]
    for at,i in cands:
        if all(abs(i-j)>20000 for _,j in kept): kept.append((at,i))
        if len(kept)>=2: break
    print(f"{ch}\t{kept[0][1]}\t{kept[0][0]:.4f}\t{kept[1][1]}\t{kept[1][0]:.4f}")
