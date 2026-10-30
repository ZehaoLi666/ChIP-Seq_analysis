fa="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
gff="/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"
seqs={};name=None;buf=[]
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name]="".join(buf)
        name=line[1:].split()[0]; buf=[]
    else: buf.append(line.strip())
seqs[name]="".join(buf)

genes={}
for line in open(gff):
    if line.startswith("#"): continue
    f=line.split("\t")
    if len(f)>4 and f[2]=="gene":
        genes.setdefault(f[0],[]).append((int(f[3])-1,int(f[4])))
for k in genes: genes[k].sort()

def gap_at(chrom, pos):
    """size of the intergenic gap containing pos, 0 if inside a gene"""
    gl=genes.get(chrom,[])
    prev_end=0
    for s,e in gl:
        if s>pos:
            return s-prev_end
        prev_end=max(prev_end,e)
    return 0

W=2000
print("chrom\tlength\tcen_start\tcen_end\tcen_len\tAT_frac\tintergenic_gap_bp\tmid_pct_of_chrom")
rows=[]
for chrom,s in sorted(seqs.items()):
    S=s.upper(); L=len(S)
    if L<50000: continue
    cands=[]
    for i in range(60000,L-W-60000,200):
        w=S[i:i+W]; cands.append(((w.count("A")+w.count("T"))/W, i))
    cands.sort(reverse=True)
    kept=[]
    for at,i in cands:
        if all(abs(i-j)>10000 for _,j in kept): kept.append((at,i))
        if len(kept)>=8: break
    # pick the top candidate that sits in a large gene desert (>=5 kb)
    pick=None
    for at,i in kept:
        g=gap_at(chrom,i+W//2)
        if g>=5000: pick=(at,i,g); break
    if pick is None: pick=(kept[0][0],kept[0][1],gap_at(chrom,kept[0][1]+W//2))
    at,i,g=pick
    a,b=i,i+W
    while a-100>0 and (S[a-100:a].count("A")+S[a-100:a].count("T"))/100>=0.92: a-=100
    while b+100<L and (S[b:b+100].count("A")+S[b:b+100].count("T"))/100>=0.92: b+=100
    fat=(S[a:b].count("A")+S[a:b].count("T"))/(b-a)
    mid=(a+b)//2
    print(f"{chrom}\t{L}\t{a}\t{b}\t{b-a}\t{fat:.4f}\t{g}\t{100*mid/L:.1f}")
    rows.append((chrom,a,b,b-a,fat,g))
with open("/rhome/zli529/lab/Aug2026_CRWN_ChIPseq/analysis/Pf3D7_v48_centromeres.bed","w") as f:
    for chrom,a,b,ln,fat,g in rows:
        f.write(f"{chrom}\t{a}\t{b}\tCEN_{chrom.split('_')[1]}\t{fat:.4f}\t.\n")
