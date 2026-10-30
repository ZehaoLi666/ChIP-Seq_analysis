#!/bin/bash
# Heterochromatin-target enrichment for CRWN-HA peaks (HPCC-side).
#
# Tests the prediction from Hollin et al. (Annu Rev Microbiol) that a
# CRWN-like nuclear-lamina protein colocalizes with H3K9me3, i.e. with
# subtelomeric regions and variant surface antigen (VSA) gene families.
#
# Emits the counts consumed by 06_enrichment_tests.py.
set -euo pipefail

PROJ=/rhome/zli529/lab/Aug2026_CRWN_ChIPseq
GD=/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
BEDTOOLS=$PKG/bedtools/2.30.0/bin/bedtools
GFF=$GD/PlasmoDB-48_Pfalciparum3D7.gff
CS=$GD/genome.chrom.sizes

cd "$PROJ"
mkdir -p analysis

# ---- 1. distance of every peak midpoint to the nearest chromosome end ----
for CTL in IGG INPUT; do
  awk 'BEGIN{OFS="\t"} NR==FNR{len[$1]=$2; next}
       ($1 in len){mid=int(($2+$3)/2); d=(mid<len[$1]-mid)?mid:len[$1]-mid;
                   print $1,$2,$3,$3-$2,d,len[$1]}' \
      "$CS" peaks/CRWN_HA_${CTL}_broad.merged_union.bed \
    > analysis/peakdist_${CTL}_broad.tsv
  awk 'BEGIN{OFS="\t"} NR==FNR{len[$1]=$2; next}
       ($1 in len){mid=int(($2+$3)/2); d=(mid<len[$1]-mid)?mid:len[$1]-mid;
                   print $1,$2,$3,$3-$2,d,len[$1]}' \
      "$CS" peaks/CRWN_HA_${CTL}_broad.reproducible_3of3.bed \
    > analysis/peakdist_${CTL}_broad_3of3.tsv
done

# expected genome fraction within X kb of either chromosome end
echo "=== EXPECTED (genome fraction near ends) ==="
for K in 50 100 200; do
  awk -v k="$K" '{w=2*k*1000; if(w>$2) w=$2; tot+=w; g+=$2}
                 END{printf "within_%dkb_genome_frac=%.4f\n", k, tot/g}' "$CS"
done

echo "=== OBSERVED (peak fraction near ends) ==="
for F in peakdist_IGG_broad peakdist_INPUT_broad \
         peakdist_IGG_broad_3of3 peakdist_INPUT_broad_3of3; do
  awk -v f="$F" '{n++; if($5<=50000)a++; if($5<=100000)b++; if($5<=200000)cc++}
                 END{printf "%s\tn=%d\tfrac50=%.4f\tfrac100=%.4f\tfrac200=%.4f\n",
                     f, n, a/n, b/n, cc/n}' "analysis/$F.tsv"
done

# ---- 2. variant surface antigen gene families ----
# NOTE: v48 descriptions label var genes "PfEMP1", NOT "var" or "erythrocyte
# membrane protein 1" -- matching on those terms silently returns zero hits.
awk 'BEGIN{OFS="\t"} $3=="gene" {
   fam="other";
   if ($0 ~ /PfEMP1/)      fam="var";
   else if ($0 ~ /rifin/)  fam="rifin";
   else if ($0 ~ /stevor/) fam="stevor";
   id="NA"; if (match($0,/ID=[^;]+/)) id=substr($0,RSTART+3,RLENGTH-3);
   print $1,$4-1,$5,id,fam,$7 }' "$GFF" \
  | sort -k1,1 -k2,2n > analysis/genes_fam.bed
awk '$5!="other"' analysis/genes_fam.bed > analysis/vsa_genes.bed

echo "=== VSA FAMILY COUNTS ==="
awk '{c[$5]++} END{for(k in c) print k, c[k]}' analysis/genes_fam.bed | sort

VSABP=$(awk '{s+=$3-$2} END{print s}' analysis/vsa_genes.bed)
GBP=$(awk '{s+=$2} END{print s}' "$CS")
awk -v v="$VSABP" -v g="$GBP" \
    'BEGIN{printf "VSA_bp=%d genome_bp=%d expected_frac=%.4f\n", v, g, v/g}'

echo "=== PEAK OVERLAP WITH VSA GENES ==="
for CTL in IGG INPUT; do
  for PK in peaks/CRWN_HA_${CTL}_broad.merged_union.bed \
            peaks/CRWN_HA_${CTL}_broad.reproducible_3of3.bed; do
    n=$(wc -l < "$PK")
    ov=$($BEDTOOLS intersect -a "$PK" -b analysis/vsa_genes.bed -u | wc -l)
    awk -v p="$PK" -v n="$n" -v o="$ov" \
        'BEGIN{printf "%s\tn=%d\tVSA_overlap=%d\tfrac=%.4f\n", p, n, o, o/n}'
  done
done

# ---- 3. centromere proximity (uses AT-derived centromere BED) ----
awk 'BEGIN{OFS="\t"}{s=$2-5000; if(s<0)s=0; print $1,s,$3+5000,$4}' \
    analysis/Pf3D7_v48_centromeres.bed > analysis/cen_flank5kb.bed

echo "=== PEAK OVERLAP WITH CENTROMERES (+/- 5 kb) ==="
for CTL in IGG INPUT; do
  for PK in peaks/CRWN_HA_${CTL}_broad.merged_union.bed \
            peaks/CRWN_HA_${CTL}_broad.reproducible_3of3.bed; do
    n=$(wc -l < "$PK")
    ov=$($BEDTOOLS intersect -a "$PK" -b analysis/cen_flank5kb.bed -u | wc -l)
    echo -e "$(basename "$PK")\tn=$n\tcen_overlap=$ov"
  done
done
