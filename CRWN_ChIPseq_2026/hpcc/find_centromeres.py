import sys
fa = "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
seqs, name, buf = {}, None, []
for line in open(fa):
    if line.startswith(">"):
        if name: seqs[name] = "".join(buf)
        name = line[1:].split()[0]; buf = []
    else:
        buf.append(line.strip())
if name: seqs[name] = "".join(buf)

W, STEP = 2000, 200
print("chrom\tlength\tAT_peak_start\tAT_peak_end\tAT_frac\tgenome_AT")
rows = []
for chrom, s in sorted(seqs.items()):
    S = s.upper(); L = len(S)
    if L < 50000:   # skip apicoplast/mito
        continue
    gAT = (S.count("A")+S.count("T"))/max(1,(L - S.count("N")))
    best = (0.0, 0)
    # exclude 60 kb at each end so subtelomeric AT-rich repeats don't win
    for i in range(60000, L-W-60000, STEP):
        w = S[i:i+W]
        at = (w.count("A")+w.count("T"))/W
        if at > best[0]:
            best = (at, i)
    at, i = best
    # extend while local AT stays >= 0.93
    a, b = i, i+W
    while a-100 > 0 and (S[a-100:a].count("A")+S[a-100:a].count("T"))/100 >= 0.93: a -= 100
    while b+100 < L and (S[b:b+100].count("A")+S[b:b+100].count("T"))/100 >= 0.93: b += 100
    finalAT = (S[a:b].count("A")+S[a:b].count("T"))/(b-a)
    print(f"{chrom}\t{L}\t{a}\t{b}\t{finalAT:.4f}\t{gAT:.4f}")
    rows.append((chrom,a,b))
with open("/rhome/zli529/lab/Aug2026_CRWN_ChIPseq/analysis/centromeres_ATderived.bed","w") as f:
    for chrom,a,b in rows:
        f.write(f"{chrom}\t{a}\t{b}\tCEN_{chrom.split('_')[1]}\n")
