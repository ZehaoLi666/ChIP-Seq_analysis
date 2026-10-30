#!/usr/bin/env python3
"""
Build the Word mapping report (methods + mapping statistics + Figure 1).

Scope is deliberately limited to read mapping and mapping statistics plus the
per-chromosome coverage figure; peak calling and QC-correlation results are
NOT included.

Usage:
  python 09_mapping_report_docx.py \
      --stats CRWN_mapping_QC_summary.csv \
      --figure CRWN_coverage_panelA.png \
      --out CRWN_ChIPseq_mapping_report.docx
"""

import argparse

import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor

GREY = RGBColor(0x55, 0x55, 0x55)

COLS = ["Sample", "Type", "Read pairs in", "Overall alignment rate (%)",
        "Concordant unique (%)", "Reads after MAPQ30 filter",
        "Duplicates removed", "Duplicate rate (%)",
        "Final reads (deduplicated)", "Estimated library size"]
KEYS = ["sample", "type", "pairs_in", "overall_align_rate_pct",
        "concordant_unique_pct", "reads_after_mapq30", "duplicates_removed",
        "dup_rate_pct", "final_reads_rmdup", "est_library_size"]

METHODS = [
    ("Reference genome",
     "Reads were aligned to the P. falciparum 3D7 genome, PlasmoDB release 48 "
     "(PlasmoDB-48_Pfalciparum3D7_Genome.fasta), using a pre-built Bowtie2 index "
     "(Pfalciparum_v48). Genome size 23,332,839 bp across 14 nuclear chromosomes "
     "plus apicoplast and mitochondrial contigs."),
    ("Quality control and trimming",
     "Raw FASTQ files were trimmed with fastp v0.23.2 in paired-end mode (-w 8), "
     "using default adapter auto-detection and quality filtering. Per-sample HTML "
     "and JSON reports were retained."),
    ("Alignment",
     "Trimmed read pairs were aligned with Bowtie2 v2.5.3 using the "
     "--very-sensitive preset and 8 threads. Alignments were streamed directly "
     "into samtools v1.19.2 for coordinate sorting and indexing."),
    ("Post-alignment filtering",
     "Aligned BAMs were filtered with samtools view to retain only properly "
     "paired, primary alignments with MAPQ \u2265 30 (-q 30 -f 2 -F 3844), "
     "excluding unmapped reads, secondary and supplementary alignments, and "
     "QC-failed reads."),
    ("Duplicate removal",
     "PCR and optical duplicates were removed with the samtools name-sort \u2192 "
     "fixmate -m \u2192 coordinate-sort \u2192 markdup -r workflow. "
     "Duplicate-removed BAMs (*.rmdup.bam) are the final alignment files used "
     "for all downstream analysis."),
    ("Coverage tracks",
     "Coverage was computed from the deduplicated BAMs with deepTools 3.5.1. "
     "Per-sample CPM-normalized tracks were generated with bamCoverage (10-bp "
     "bins), and IP-minus-control tracks with bamCompare (--operation subtract, "
     "--normalizeUsing CPM, --scaleFactorsMethod None). For Figure 1 the "
     "CRWN-HA D1 and D2 IP-minus-Input tracks were averaged at 10-bp resolution "
     "and negative values were then set to zero; the averaged signal was "
     "summarized in 5-kb bins for display."),
    ("Compute",
     "Trimming and alignment were run as a five-task SLURM array (partition "
     "epyc, 8 CPUs, 60 GB RAM, 48 h wall limit; job 27547872). All five tasks "
     "completed successfully (exit code 0:0), with per-sample runtimes of "
     "1 h 08 m to 3 h 08 m."),
]

NOTES = [
    "Alignment rates for the three CRWN-HA replicates were consistent "
    "(88.25\u201389.16%). Input aligned highest (95.19%) and IgG lowest (80.47%).",
    "Duplicate rates were high across all libraries (36.8\u201353.0%). Estimated "
    "library sizes (15.9\u201334.4 M) indicate the libraries are "
    "complexity-limited rather than under-sequenced; duplicates were removed "
    "following the standard protocol.",
    "Final deduplicated depths were 44.1\u201352.4 M reads for the IP replicates, "
    "26.4 M for IgG and 28.5 M for Input.",
    "Centromere positions in Figure 1 were derived from sequence AT content, not "
    "from the PlasmoDB v48 annotation, which contains no centromere features. "
    "The chromosome 10 position is provisional: the coverage maximum lies "
    "~350 kb from the inferred coordinate.",
]

FIG_CAPTION = (
    "Figure 1. Averaged CRWN-HA read coverage along each of the 14 nuclear "
    "chromosomes. Signal is the mean of the CRWN-HA D1 and D2 IP-minus-Input "
    "tracks with negative values set to zero, summarized in 5-kb bins. All "
    "chromosomes share the same horizontal scale (absolute position in Mb) and "
    "the same vertical scale (0\u20135.8 CPM, shown on chr01); chromosome lengths "
    "are annotated at the right of each track. Red triangles mark centromere "
    "positions inferred from local AT content.")


def small_grey(paragraph, text, size=8):
    run = paragraph.add_run(text)
    run.font.size = Pt(size)
    run.font.color.rgb = GREY
    return run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--figure", required=True)
    ap.add_argument("--out", default="CRWN_ChIPseq_mapping_report.docx")
    args = ap.parse_args()

    df = pd.read_csv(args.stats)

    doc = Document()
    for section in doc.sections:
        section.left_margin = section.right_margin = Inches(0.9)
    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(10.5)

    doc.add_heading("CRWN ChIP-seq: Read Mapping and Mapping Statistics", level=1)
    sub = doc.add_paragraph()
    run = sub.add_run("Plasmodium falciparum 3D7 \u00b7 PlasmoDB v48 \u00b7 UCR HPCC")
    run.italic = True
    run.font.size = Pt(9.5)
    run.font.color.rgb = GREY

    doc.add_heading("1. Samples", level=2)
    doc.add_paragraph(
        "Five paired-end libraries sequenced on lane L008: three CRWN-HA "
        "immunoprecipitation replicates (D1, D2, D3), one IgG mock-IP control, "
        "and one Input chromatin control.")

    doc.add_heading("2. Methods", level=2)
    for heading, body in METHODS:
        doc.add_paragraph(heading, style="Heading 3")
        doc.add_paragraph(body)

    doc.add_heading("3. Mapping statistics", level=2)
    table = doc.add_table(rows=1, cols=len(COLS))
    table.style = "Light Grid Accent 1"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, name in enumerate(COLS):
        cell = table.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(name)
        run.bold = True
        run.font.size = Pt(8)
    for _, row in df.iterrows():
        cells = table.add_row().cells
        for i, key in enumerate(KEYS):
            value = row[key]
            text = (f"{value:,}" if isinstance(value, (int,))
                    or (hasattr(value, "dtype") and "int" in str(value.dtype))
                    else str(value))
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(text)
            run.font.size = Pt(8)
            if i > 1:
                cells[i].paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.RIGHT
    small_grey(doc.add_paragraph(),
               "Read pairs in = pairs surviving fastp trimming (Bowtie2 input). "
               "Overall alignment rate and concordant unique are reported by "
               "Bowtie2. Reads after MAPQ30 filter, duplicates removed, and final "
               "reads are single-read counts from samtools; duplicate rate = "
               "duplicates removed / reads after MAPQ30 filter. Estimated library "
               "size is the samtools markdup complexity estimate.")

    doc.add_heading("4. Genome-wide read coverage", level=2)
    pic = doc.add_paragraph()
    pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pic.add_run().add_picture(args.figure, width=Inches(6.4))
    small_grey(doc.add_paragraph(), FIG_CAPTION)

    doc.add_heading("5. Notes", level=2)
    for note in NOTES:
        doc.add_paragraph(note, style="List Bullet")

    doc.save(args.out)

    check = Document(args.out)
    print(f"saved {args.out}: paragraphs={len(check.paragraphs)} "
          f"tables={len(check.tables)} images={len(check.inline_shapes)}")


if __name__ == "__main__":
    main()
