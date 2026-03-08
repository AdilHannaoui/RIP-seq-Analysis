# =============================================================================
# Snakefile — RIP-seq Analysis Pipeline
# Author: Adil Hannaoui Anaaoui
#
# Pipeline:
#   FASTQ → FastQC (pre) → Trimmomatic + FastQC (post) → Bowtie2 + BedGraph
#       → featureCounts
#       → MACS3 peak calling (per IP/IN pair)
#       → Peak intersection across replicates (per condition)
#       → TSS window analysis (merged IP BAMs)
#       → TSS metagene (bigWig + deepTools)
#
# Triple-evidence strategy:
#   Genes appearing in featureCounts + MACS3 + window analysis
#   are considered high-confidence Rpb4-bound transcripts.
#
# Usage:
#   snakemake --cores 10
#   snakemake --dag | dot -Tsvg > dag.svg
# =============================================================================

from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
configfile: "config/config.yaml"

WORKDIR      = config["workdir"]
FASTQ_DIR    = config["fastq_dir"]
OUTPUT_DIR   = config["output_dir"]
PLOTS_DIR    = config["plots_dir"]
BOWTIE2_IDX  = config["bowtie2_index"]
GTF          = config["gtf"]
GENOME_FA    = config["genome_fa"]
TRIMMO_JAR   = config["trimmomatic_jar"]
ADAPTERS     = config["adapters"]
THREADS      = config["threads"]

# ── Samples ───────────────────────────────────────────────────────────────────
CONDITIONS = ["WT", "M1", "M12"]
TYPES      = ["IP", "IN"]
REPS       = ["1", "2", "3"]

# All 18 samples: WT_IP1, WT_IN1, M1_IP1, ...
SAMPLES = [
    f"{cond}_{t}{rep}"
    for cond in CONDITIONS
    for t in TYPES
    for rep in REPS
]

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        # FastQC pre
        expand(f"{OUTPUT_DIR}/fastqc_pre/{{sample}}_fastqc.html", sample=SAMPLES),
        # FastQC post
        expand(f"{OUTPUT_DIR}/fastqc_post/{{sample}}_trimmed_fastqc.html", sample=SAMPLES),
        # featureCounts
        expand(f"{OUTPUT_DIR}/featurecounts/Counts_{{sample}}.txt", sample=SAMPLES),
        # MACS3 peak calling (one per IP/IN pair)
        expand(
            f"{OUTPUT_DIR}/macs3/{{condition}}_rep{{rep}}_peaks.narrowPeak",
            condition=CONDITIONS, rep=REPS
        ),
        # Peak intersection (one per condition)
        expand(
            f"{OUTPUT_DIR}/macs3/{{condition}}_common_counts.txt",
            condition=CONDITIONS
        ),
        # Window analysis (one per condition, merged IP)
        expand(
            f"{OUTPUT_DIR}/Windows_Analysis/{{condition}}_IP_merged_sorted_TSS_1kb_coverage.bed",
            condition=CONDITIONS
        ),
        # TSS metagene plots
        f"{PLOTS_DIR}/TSS_metagene_IP.png",
        f"{PLOTS_DIR}/TSS_metagene_IN.png",


# ── Rule 1: FastQC pre-trimming ───────────────────────────────────────────────
rule fastqc_pre:
    input:
        fastq = f"{FASTQ_DIR}/{{sample}}.fastq"
    output:
        html = f"{OUTPUT_DIR}/fastqc_pre/{{sample}}_fastqc.html",
        zip  = f"{OUTPUT_DIR}/fastqc_pre/{{sample}}_fastqc.zip"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_pre.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/fastqc_pre
        fastqc {input.fastq} --threads {threads} --outdir {OUTPUT_DIR}/fastqc_pre > {log} 2>&1
        """


# ── Rule 2: Trimmomatic ───────────────────────────────────────────────────────
rule trimming:
    input:
        fastq = f"{FASTQ_DIR}/{{sample}}.fastq"
    output:
        trimmed = f"{OUTPUT_DIR}/fastq_trimmed/{{sample}}_trimmed.fastq"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_trimming.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/fastq_trimmed
        java -jar {TRIMMO_JAR} SE -threads {threads} \
            {input.fastq} {output.trimmed} \
            ILLUMINACLIP:{ADAPTERS}:2:30:10 \
            SLIDINGWINDOW:4:20 MINLEN:20 -phred33 \
            > {log} 2>&1
        """


# ── Rule 3: FastQC post-trimming ──────────────────────────────────────────────
rule fastqc_post:
    input:
        trimmed = f"{OUTPUT_DIR}/fastq_trimmed/{{sample}}_trimmed.fastq"
    output:
        html = f"{OUTPUT_DIR}/fastqc_post/{{sample}}_trimmed_fastqc.html",
        zip  = f"{OUTPUT_DIR}/fastqc_post/{{sample}}_trimmed_fastqc.zip"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_fastqc_post.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/fastqc_post
        fastqc {input.trimmed} --threads {threads} --outdir {OUTPUT_DIR}/fastqc_post > {log} 2>&1
        """


# ── Rule 4: Bowtie2 alignment + BAM index + BedGraph ─────────────────────────
rule bowtie2:
    input:
        trimmed = f"{OUTPUT_DIR}/fastq_trimmed/{{sample}}_trimmed.fastq"
    output:
        bam      = f"{OUTPUT_DIR}/bowtie2/{{sample}}.bam",
        bai      = f"{OUTPUT_DIR}/bowtie2/{{sample}}.bam.bai",
        bedgraph = f"{OUTPUT_DIR}/visualization/{{sample}}.bedgraph"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_bowtie2.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/bowtie2 {OUTPUT_DIR}/visualization
        bowtie2 -p {threads} --very-fast \
            -x {BOWTIE2_IDX} \
            -U {input.trimmed} \
            2> {log} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}

        bedtools genomecov -ibam {output.bam} -bg \
        | sort -k1,1 -k2,2n > {output.bedgraph}
        """


# ── Rule 5: featureCounts ─────────────────────────────────────────────────────
rule featurecounts:
    input:
        bam = f"{OUTPUT_DIR}/bowtie2/{{sample}}.bam"
    output:
        full   = f"{OUTPUT_DIR}/featurecounts/{{sample}}_featurecounts.txt",
        counts = f"{OUTPUT_DIR}/featurecounts/Counts_{{sample}}.txt"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}_featurecounts.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/featurecounts
        featureCounts -T {threads} -s 2 \
            -a {GTF} -o {output.full} {input.bam} \
            > {log} 2>&1
        awk 'NR > 1 && !/^#/ {{print $1"\t"$7}}' {output.full} > {output.counts}
        """


# ── Rule 6: MACS3 peak calling (per condition + replicate) ───────────────────
rule macs3:
    input:
        ip = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IP{{rep}}.bam",
        in_ = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IN{{rep}}.bam"
    output:
        peaks = f"{OUTPUT_DIR}/macs3/{{condition}}_rep{{rep}}_peaks.narrowPeak"
    log:
        f"{OUTPUT_DIR}/logs/macs3_{{condition}}_rep{{rep}}.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/macs3
        macs3 callpeak \
            -t {input.ip} \
            -c {input.in_} \
            --format BAM \
            --name {wildcards.condition}_rep{wildcards.rep} \
            --outdir {OUTPUT_DIR}/macs3 \
            --pvalue 0.1 \
            --call-summits \
            --nomodel \
            --extsize 150 \
            --gsize 1.2e7 \
            --keep-dup all \
            > {log} 2>&1
        """


# ── Rule 7: Peak intersection across replicates (per condition) ───────────────
rule peak_intersection:
    input:
        peaks = expand(
            f"{OUTPUT_DIR}/macs3/{{{{condition}}}}_rep{{rep}}_peaks.narrowPeak",
            rep=REPS
        ),
        bams = expand(
            f"{OUTPUT_DIR}/bowtie2/{{{{condition}}}}_{t}{{rep}}.bam",
            t=TYPES, rep=REPS
        )
    output:
        counts = f"{OUTPUT_DIR}/macs3/{{condition}}_common_counts.txt"
    log:
        f"{OUTPUT_DIR}/logs/{{condition}}_intersection.log"
    threads: 1
    shell:
        """
        # Chain intersect: rep1 ∩ rep2 ∩ rep3
        cp {input.peaks[0]} {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed

        bedtools intersect \
            -a {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed \
            -b {input.peaks[1]} -wa -u \
            > {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.tmp
        mv {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.tmp \
           {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed

        bedtools intersect \
            -a {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed \
            -b {input.peaks[2]} -wa -u \
            > {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.tmp
        mv {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.tmp \
           {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed

        # Annotate with GTF
        bedtools intersect \
            -a {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed \
            -b {GTF} -wa -wb \
        | awk -F'\\t' '{{
            gene_id="NA"; biotype="NA"
            if (match($0, /gene_id "([^"]+)"/, m)) gene_id=m[1]
            if (match($0, /gene_biotype "([^"]+)"/, b)) biotype=b[1]
            print $1"\\t"$2"\\t"$3"\\t."\\t"."\\t"."\\t"gene_id"\\t"biotype
        }}' | sort -u \
        > {OUTPUT_DIR}/macs3/{wildcards.condition}_annotated.bed

        # multicov: count reads per IP/IN BAM across annotated peaks
        echo -e "chrom\\tstart\\tend\\tname\\tscore\\tstrand\\tGene_id\\tBiotype\\tIP1\\tIN1\\tIP2\\tIN2\\tIP3\\tIN3" \
            > {output.counts}

        bedtools multicov \
            -bams {input.bams} \
            -bed {OUTPUT_DIR}/macs3/{wildcards.condition}_annotated.bed \
            >> {output.counts}

        # Cleanup intermediates
        rm -f {OUTPUT_DIR}/macs3/{wildcards.condition}_raw.bed \
              {OUTPUT_DIR}/macs3/{wildcards.condition}_annotated.bed

        echo "Common peaks for {wildcards.condition}: $(tail -n +2 {output.counts} | wc -l)" >> {log}
        """


# ── Rule 8: Merge IP replicates per condition ─────────────────────────────────
rule merge_ip_bams:
    input:
        bams = expand(
            f"{OUTPUT_DIR}/bowtie2/{{{{condition}}}}_IP{{rep}}.bam",
            rep=REPS
        )
    output:
        merged = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IP_merged_sorted.bam",
        bai    = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IP_merged_sorted.bam.bai"
    log:
        f"{OUTPUT_DIR}/logs/{{condition}}_merge_ip.log"
    threads: THREADS
    shell:
        """
        samtools merge -f {OUTPUT_DIR}/bowtie2/{wildcards.condition}_IP_merged.bam {input.bams}
        samtools sort -@ {threads} \
            -o {output.merged} \
            {OUTPUT_DIR}/bowtie2/{wildcards.condition}_IP_merged.bam
        samtools index {output.merged}
        rm -f {OUTPUT_DIR}/bowtie2/{wildcards.condition}_IP_merged.bam
        """


# ── Rule 9: TSS window analysis ───────────────────────────────────────────────
rule tss_windows:
    input:
        merged = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IP_merged_sorted.bam",
        bai    = f"{OUTPUT_DIR}/bowtie2/{{condition}}_IP_merged_sorted.bam.bai"
    output:
        coverage = f"{OUTPUT_DIR}/Windows_Analysis/{{condition}}_IP_merged_sorted_TSS_1kb_coverage.bed"
    log:
        f"{OUTPUT_DIR}/logs/{{condition}}_window_analysis.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTPUT_DIR}/Windows_Analysis

        # Extract TSS from GTF (strand-aware)
        awk '$3=="gene" {{
            if ($7=="+") {{ tss=$4-1 }} else {{ tss=$5-1 }}
            if (tss>=0) {{
                match($0, /gene_id "([^"]+)"/, m)
                print $1"\\t"tss"\\t"tss+1"\\t"m[1]"\\t0\\t"$7
            }}
        }}' {GTF} > {OUTPUT_DIR}/Windows_Analysis/genes_TSS.bed

        # Genome sizes from FASTA index
        if [[ ! -f {GENOME_FA}.fai ]]; then
            samtools faidx {GENOME_FA}
        fi
        cut -f1,2 {GENOME_FA}.fai > {OUTPUT_DIR}/Windows_Analysis/genome.sizes

        # ±500 bp windows around TSS
        bedtools slop \
            -b 500 \
            -g {OUTPUT_DIR}/Windows_Analysis/genome.sizes \
            -i {OUTPUT_DIR}/Windows_Analysis/genes_TSS.bed \
            > {OUTPUT_DIR}/Windows_Analysis/genes_TSS_1kb.bed

        # Coverage in windows
        bedtools coverage \
            -a {OUTPUT_DIR}/Windows_Analysis/genes_TSS_1kb.bed \
            -b {input.merged} \
            -counts \
            > {output.coverage}
        """


# ── Rule 10: TSS metagene (bigWig + deepTools) ────────────────────────────────
rule tss_metagene:
    input:
        bedgraphs_ip = expand(
            f"{OUTPUT_DIR}/visualization/{{condition}}_IP{{rep}}.bedgraph",
            condition=CONDITIONS, rep=REPS
        ),
        bedgraphs_in = expand(
            f"{OUTPUT_DIR}/visualization/{{condition}}_IN{{rep}}.bedgraph",
            condition=CONDITIONS, rep=REPS
        )
    output:
        plot_ip = f"{PLOTS_DIR}/TSS_metagene_IP.png",
        plot_in = f"{PLOTS_DIR}/TSS_metagene_IN.png"
    log:
        f"{OUTPUT_DIR}/logs/tss_metagene.log"
    threads: THREADS
    shell:
        """
        mkdir -p {PLOTS_DIR} {PLOTS_DIR}/merged_bw

        # Fix chromosome names (Roman numerals → chr prefix)
        fix_bedgraph() {{
            local f="$1"
            local out="${{f%.bedgraph}}.fixed.bedgraph"
            sed -E \
                -e 's/^I\t/chrI\t/' -e 's/^II\t/chrII\t/' \
                -e 's/^III\t/chrIII\t/' -e 's/^IV\t/chrIV\t/' \
                -e 's/^V\t/chrV\t/' -e 's/^VI\t/chrVI\t/' \
                -e 's/^VII\t/chrVII\t/' -e 's/^VIII\t/chrVIII\t/' \
                -e 's/^IX\t/chrIX\t/' -e 's/^X\t/chrX\t/' \
                -e 's/^XI\t/chrXI\t/' -e 's/^XII\t/chrXII\t/' \
                -e 's/^XIII\t/chrXIII\t/' -e 's/^XIV\t/chrXIV\t/' \
                -e 's/^XV\t/chrXV\t/' -e 's/^XVI\t/chrXVI\t/' \
                -e 's/^Mito\t/chrM\t/' \
                "$f" | awk 'NF==4' | sort -k1,1 -k2,2n > "$out"
            echo "$out"
        }}
        export -f fix_bedgraph

        # Chromosome sizes from fixed bedGraphs
        ALL_FIXED=()
        for f in {input.bedgraphs_ip} {input.bedgraphs_in}; do
            fixed=$(fix_bedgraph "$f")
            ALL_FIXED+=("$fixed")
        done

        CHROM_SIZES="{PLOTS_DIR}/chrom.sizes"
        cat "${{ALL_FIXED[@]}}" \
        | awk '{{if($3>max[$1]) max[$1]=$3}} END{{for(c in max) print c"\\t"max[c]}}' \
        | sort -k1,1 > "$CHROM_SIZES"

        # bedGraph → bigWig via pyBigWig
        bg2bw() {{
            python - "$1" "$CHROM_SIZES" "${{1%.fixed.bedgraph}}.bw" << 'PYEOF'
import pyBigWig, sys
bg, cs, bw_out = sys.argv[1], sys.argv[2], sys.argv[3]
sizes = {{l.split()[0]: int(l.split()[1]) for l in open(cs)}}
bw = pyBigWig.open(bw_out, "w")
bw.addHeader(list(sizes.items()))
for line in open(bg):
    if line.startswith(('#','track')): continue
    c,s,e,v = line.strip().split()
    bw.addEntries([c],[int(s)],ends=[int(e)],values=[float(v)])
bw.close()
PYEOF
        }}
        export -f bg2bw

        for f in "${{ALL_FIXED[@]}}"; do bg2bw "$f"; done

        # Merge replicates per condition using bigwigCompare (mean)
        for cond in {CONDITIONS}; do
            for type in IP IN; do
                bigwigCompare \
                    -b1 {OUTPUT_DIR}/visualization/${{cond}}_${{type}}1.bw \
                    -b2 {OUTPUT_DIR}/visualization/${{cond}}_${{type}}2.bw \
                    --operation mean -o /tmp/${{cond}}_${{type}}_temp.bw -p {threads}
                bigwigCompare \
                    -b1 /tmp/${{cond}}_${{type}}_temp.bw \
                    -b2 {OUTPUT_DIR}/visualization/${{cond}}_${{type}}3.bw \
                    --operation mean \
                    -o {PLOTS_DIR}/merged_bw/${{cond}}_${{type}}_merged.bw \
                    -p {threads}
            done
        done
        rm -f /tmp/*_temp.bw

        # Fix GTF chromosome names
        GTF_FIXED="{OUTPUT_DIR}/gtf_fixed.gtf"
        sed -E \
            -e 's/^I\t/chrI\t/' -e 's/^II\t/chrII\t/' \
            -e 's/^III\t/chrIII\t/' -e 's/^IV\t/chrIV\t/' \
            -e 's/^V\t/chrV\t/' -e 's/^VI\t/chrVI\t/' \
            -e 's/^Mito\t/chrM\t/' {GTF} > "$GTF_FIXED"

        # TSS BED from fixed GTF
        TSS_BED="{PLOTS_DIR}/TSS_regions.bed"
        awk '$3=="gene" {{
            match($0, /gene_id "([^"]+)"/, m)
            if ($7=="+") {{ s=$4-1; e=$4 }} else {{ s=$5-1; e=$5 }}
            if (s<0) s=0
            print $1"\\t"s"\\t"e"\\t"m[1]"\\t0\\t"$7
        }}' "$GTF_FIXED" > "$TSS_BED"

        # IP merged bigWigs
        IP_BWS=({PLOTS_DIR}/merged_bw/WT_IP_merged.bw \
                {PLOTS_DIR}/merged_bw/M1_IP_merged.bw \
                {PLOTS_DIR}/merged_bw/M12_IP_merged.bw)

        IN_BWS=({PLOTS_DIR}/merged_bw/WT_IN_merged.bw \
                {PLOTS_DIR}/merged_bw/M1_IN_merged.bw \
                {PLOTS_DIR}/merged_bw/M12_IN_merged.bw)

        # computeMatrix + plotProfile IP
        computeMatrix reference-point \
            --referencePoint TSS -b 1000 -a 1000 \
            -R "$TSS_BED" -S "${{IP_BWS[@]}}" \
            --binSize 10 --skipZeros \
            -o {PLOTS_DIR}/TSS_matrix_IP.gz -p {threads}

        plotProfile \
            -m {PLOTS_DIR}/TSS_matrix_IP.gz \
            -o {output.plot_ip} \
            --plotTitle "TSS Metagene - IP Samples" \
            --yAxisLabel "Normalized signal" \
            --colors "#E41A1C" "#377EB8" "#4DAF4A" \
            --plotHeight 10 --plotWidth 15 --dpi 300

        # computeMatrix + plotProfile IN
        computeMatrix reference-point \
            --referencePoint TSS -b 1000 -a 1000 \
            -R "$TSS_BED" -S "${{IN_BWS[@]}}" \
            --binSize 10 --skipZeros \
            -o {PLOTS_DIR}/TSS_matrix_IN.gz -p {threads}

        plotProfile \
            -m {PLOTS_DIR}/TSS_matrix_IN.gz \
            -o {output.plot_in} \
            --plotTitle "TSS Metagene - IN Samples" \
            --yAxisLabel "Normalized signal" \
            --colors "#E41A1C" "#377EB8" "#4DAF4A" \
            --plotHeight 10 --plotWidth 15 --dpi 300

        # Cleanup
        rm -f "${{ALL_FIXED[@]}}" {OUTPUT_DIR}/visualization/*.bw \
              "$GTF_FIXED" "$CHROM_SIZES"
        """
