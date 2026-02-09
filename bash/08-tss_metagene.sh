#!/usr/bin/env bash
set -euo pipefail

echo "=== TSS Metagene (bedGraph to bigWig) ==="

# Cargar config
source "$(dirname "$0")/config.sh"

# Directorios
VIS_DIR="$OUTPUT_DIR/visualization"
OUT_DIR="$PLOTS_DIR"
mkdir -p "$OUT_DIR"

echo "Using bedGraphs from: $VIS_DIR"
echo "Saving results to: $OUT_DIR"

# ================================
# 0. Renombrar cromosomas en bedGraphs
# ================================
echo "Fixing chromosome names in bedGraphs..."

mapfile -t bedgraphs < <(find "$VIS_DIR" -maxdepth 1 -name "*.bedgraph" | sort)

for f in "${bedgraphs[@]}"; do
    fixed="${f%.bedgraph}.fixed.bedgraph"

    sed -E \
        -e 's/^I[[:space:]]/chrI\t/' \
        -e 's/^II[[:space:]]/chrII\t/' \
        -e 's/^III[[:space:]]/chrIII\t/' \
        -e 's/^IV[[:space:]]/chrIV\t/' \
        -e 's/^V[[:space:]]/chrV\t/' \
        -e 's/^VI[[:space:]]/chrVI\t/' \
        -e 's/^VII[[:space:]]/chrVII\t/' \
        -e 's/^VIII[[:space:]]/chrVIII\t/' \
        -e 's/^IX[[:space:]]/chrIX\t/' \
        -e 's/^X[[:space:]]/chrX\t/' \
        -e 's/^XI[[:space:]]/chrXI\t/' \
        -e 's/^XII[[:space:]]/chrXII\t/' \
        -e 's/^XIII[[:space:]]/chrXIII\t/' \
        -e 's/^XIV[[:space:]]/chrXIV\t/' \
        -e 's/^XV[[:space:]]/chrXV\t/' \
        -e 's/^XVI[[:space:]]/chrXVI\t/' \
        -e 's/^Mito/chrM/' \
        "$f" \
    | awk 'NF==4' \
    | sort -k1,1 -k2,2n \
    > "$fixed"
done

# ================================
# 1. Renombrar cromosomas en el GTF
# ================================
echo "Fixing chromosome names in GTF..."

GTF_FIXED="${GTF%.gtf}.fixed.gtf"

sed -E \
    -e 's/^I\t/chrI\t/' \
    -e 's/^II\t/chrII\t/' \
    -e 's/^III\t/chrIII\t/' \
    -e 's/^IV\t/chrIV\t/' \
    -e 's/^V\t/chrV\t/' \
    -e 's/^VI\t/chrVI\t/' \
    -e 's/^VII\t/chrVII\t/' \
    -e 's/^VIII\t/chrVIII\t/' \
    -e 's/^IX\t/chrIX\t/' \
    -e 's/^X\t/chrX\t/' \
    -e 's/^XI\t/chrXI\t/' \
    -e 's/^XII\t/chrXII\t/' \
    -e 's/^XIII\t/chrXIII\t/' \
    -e 's/^XIV\t/chrXIV\t/' \
    -e 's/^XV\t/chrXV\t/' \
    -e 's/^XVI\t/chrXVI\t/' \
    -e 's/^Mito\t/chrM\t/' \
    "$GTF" > "$GTF_FIXED"

# ================================
# 2. Detectar bedGraphs IP e IN
# ================================
echo "Detecting bedGraph files..."

mapfile -t ip_beds < <(find "$VIS_DIR" -maxdepth 1 -name "*_IP*.fixed.bedgraph" | sort)
mapfile -t in_beds < <(find "$VIS_DIR" -maxdepth 1 -name "*_IN*.fixed.bedgraph" | sort)

if (( ${#ip_beds[@]} == 0 )); then
    echo "ERROR: No IP bedGraph files found."
    exit 1
fi

if (( ${#in_beds[@]} == 0 )); then
    echo "ERROR: No IN bedGraph files found."
    exit 1
fi

echo "Found IP bedGraphs:"
printf "  %s\n" "${ip_beds[@]}"

echo "Found IN bedGraphs:"
printf "  %s\n" "${in_beds[@]}"

# ================================
# 3. Generar archivo de tamaños de cromosomas desde bedGraphs
# ================================
echo "Generating chromosome sizes from bedGraphs..."

CHROM_SIZES="$OUT_DIR/chrom.sizes"

# Combinar todos los bedGraphs y obtener el máximo end para cada cromosoma
cat "${ip_beds[@]}" "${in_beds[@]}" | \
    awk '{print $1, $3}' | \
    sort -k1,1 -k2,2n | \
    awk '{if($2>max[$1]) max[$1]=$2} END{for(chr in max) print chr"\t"max[chr]}' | \
    sort -k1,1 > "$CHROM_SIZES"

echo "Chromosome sizes saved at: $CHROM_SIZES"
echo "Content preview:"
head "$CHROM_SIZES"

# ================================
# 4. Convertir bedGraphs a bigWig (usando pyBigWig)
# ================================
echo "Converting bedGraphs to bigWig format using pyBigWig..."

# Crear script Python temporal
cat > /tmp/bg2bw.py << 'PYEOF'
import pyBigWig
import sys

bedgraph_file = sys.argv[1]
chrom_sizes_file = sys.argv[2]
output_file = sys.argv[3]

# Leer tamaños de cromosomas
chrom_sizes = {}
with open(chrom_sizes_file) as f:
    for line in f:
        chrom, size = line.strip().split()
        chrom_sizes[chrom] = int(size)

# Crear bigWig
bw = pyBigWig.open(output_file, "w")
bw.addHeader(list(chrom_sizes.items()))

# Leer bedGraph y añadir al bigWig
with open(bedgraph_file) as f:
    for line in f:
        if line.startswith('#') or line.startswith('track'):
            continue
        parts = line.strip().split()
        if len(parts) == 4:
            chrom, start, end, value = parts
            bw.addEntries([chrom], [int(start)], ends=[int(end)], values=[float(value)])

bw.close()
print(f"Created: {output_file}")
PYEOF

# Convertir IP bedGraphs
ip_bws=()
for bed in "${ip_beds[@]}"; do
    bw="${bed%.fixed.bedgraph}.bw"
    echo "  Converting: $(basename "$bed") -> $(basename "$bw")"
    python /tmp/bg2bw.py "$bed" "$CHROM_SIZES" "$bw"
    ip_bws+=("$bw")
done

# Convertir IN bedGraphs
in_bws=()
for bed in "${in_beds[@]}"; do
    bw="${bed%.fixed.bedgraph}.bw"
    echo "  Converting: $(basename "$bed") -> $(basename "$bw")"
    python /tmp/bg2bw.py "$bed" "$CHROM_SIZES" "$bw"
    in_bws+=("$bw")
done

# Limpiar script temporal
rm /tmp/bg2bw.py

echo "BigWig conversion complete."
echo "IP bigWigs:"
printf "  %s\n" "${ip_bws[@]}"
echo "IN bigWigs:"
printf "  %s\n" "${in_bws[@]}"


# ================================
# 4b. Mergear réplicas por condición
# ================================
echo "Merging replicates by condition..."

MERGED_DIR="$OUT_DIR/merged_bw"
mkdir -p "$MERGED_DIR"

# WT IP
echo "  Merging WT IP..."
bigwigCompare -b1 "$VIS_DIR/WT_IP1.bw" -b2 "$VIS_DIR/WT_IP2.bw" --operation mean -o /tmp/WT_IP_temp.bw
bigwigCompare -b1 /tmp/WT_IP_temp.bw -b2 "$VIS_DIR/WT_IP3.bw" --operation mean -o "$MERGED_DIR/WT_IP_merged.bw"

# WT IN
echo "  Merging WT IN..."
bigwigCompare -b1 "$VIS_DIR/WT_IN1.bw" -b2 "$VIS_DIR/WT_IN2.bw" --operation mean -o /tmp/WT_IN_temp.bw
bigwigCompare -b1 /tmp/WT_IN_temp.bw -b2 "$VIS_DIR/WT_IN3.bw" --operation mean -o "$MERGED_DIR/WT_IN_merged.bw"

# M1 IP
echo "  Merging M1 IP..."
bigwigCompare -b1 "$VIS_DIR/M1_IP1.bw" -b2 "$VIS_DIR/M1_IP2.bw" --operation mean -o /tmp/M1_IP_temp.bw
bigwigCompare -b1 /tmp/M1_IP_temp.bw -b2 "$VIS_DIR/M1_IP3.bw" --operation mean -o "$MERGED_DIR/M1_IP_merged.bw"

# M1 IN
echo "  Merging M1 IN..."
bigwigCompare -b1 "$VIS_DIR/M1_IN1.bw" -b2 "$VIS_DIR/M1_IN2.bw" --operation mean -o /tmp/M1_IN_temp.bw
bigwigCompare -b1 /tmp/M1_IN_temp.bw -b2 "$VIS_DIR/M1_IN3.bw" --operation mean -o "$MERGED_DIR/M1_IN_merged.bw"

# M12 IP
echo "  Merging M12 IP..."
bigwigCompare -b1 "$VIS_DIR/M12_IP1.bw" -b2 "$VIS_DIR/M12_IP2.bw" --operation mean -o /tmp/M12_IP_temp.bw
bigwigCompare -b1 /tmp/M12_IP_temp.bw -b2 "$VIS_DIR/M12_IP3.bw" --operation mean -o "$MERGED_DIR/M12_IP_merged.bw"

# M12 IN
echo "  Merging M12 IN..."
bigwigCompare -b1 "$VIS_DIR/M12_IN1.bw" -b2 "$VIS_DIR/M12_IN2.bw" --operation mean -o /tmp/M12_IN_temp.bw
bigwigCompare -b1 /tmp/M12_IN_temp.bw -b2 "$VIS_DIR/M12_IN3.bw" --operation mean -o "$MERGED_DIR/M12_IN_merged.bw"

# Limpiar temporales
rm -f /tmp/*_temp.bw

# Actualizar arrays
mapfile -t ip_bws < <(find "$MERGED_DIR" -name "*_IP_merged.bw" | sort)
mapfile -t in_bws < <(find "$MERGED_DIR" -name "*_IN_merged.bw" | sort)

echo "Merged bigWigs created:"
echo "IP samples:"
printf "  %s\n" "${ip_bws[@]}"
echo "IN samples:"
printf "  %s\n" "${in_bws[@]}"

# ================================
# 5. Generar BED de TSS desde GTF
# ================================
echo "Extracting TSS from GTF..."

TSS_BED="$OUT_DIR/TSS_regions.bed"

awk '$3=="gene" {
    split($9,a,";");
    for(i in a){
        if(a[i] ~ /gene_id/){
            gsub(/gene_id |"|;/,"",a[i]);
            gene=a[i]
        }
    }
    if($7=="+"){
        start=$4-1;
        end=$4;
    } else {
        start=$5-1;
        end=$5;
    }
    if(start<0) start=0;
    print $1"\t"start"\t"end"\t"gene"\t0\t"$7
}' "$GTF_FIXED" > "$TSS_BED"

echo "TSS BED saved at: $TSS_BED"
echo "Total TSS regions: $(wc -l < "$TSS_BED")"

# ================================
# 6. computeMatrix - Crear matrices separadas
# ================================
echo "Running computeMatrix..."

MATRIX_IP="$OUT_DIR/TSS_matrix_IP.gz"
MATRIX_IN="$OUT_DIR/TSS_matrix_IN.gz"

# Matriz solo para IP
echo "  Creating IP matrix..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "$TSS_BED" \
    -S "${ip_bws[@]}" \
    --binSize 10 \
    --skipZeros \
    -o "$MATRIX_IP" \
    -p max

echo "  IP Matrix saved at: $MATRIX_IP"

# Matriz solo para IN
echo "  Creating IN matrix..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R "$TSS_BED" \
    -S "${in_bws[@]}" \
    --binSize 10 \
    --skipZeros \
    -o "$MATRIX_IN" \
    -p max

echo "  IN Matrix saved at: $MATRIX_IN"

# ================================
# 7. plotProfile - Generar gráficas separadas
# ================================
echo "Generating TSS metagene plots..."

# Gráfica para IP samples
echo "  Creating IP metagene plot..."
plotProfile \
    -m "$MATRIX_IP" \
    -o "$OUT_DIR/TSS_metagene_IP.png" \
    --plotTitle "TSS Metagene - IP Samples" \
    --regionsLabel "TSS" \
    --yAxisLabel "Normalized signal" \
    --plotHeight 10 \
    --plotWidth 15 \
    --colors "#E41A1C" "#377EB8" "#4DAF4A" \
    --plotFileFormat png \
    --dpi 300

# Gráfica para IN samples
echo "  Creating IN metagene plot..."
plotProfile \
    -m "$MATRIX_IN" \
    -o "$OUT_DIR/TSS_metagene_IN.png" \
    --plotTitle "TSS Metagene - IN Samples" \
    --regionsLabel "TSS" \
    --yAxisLabel "Normalized signal" \
    --plotHeight 10 \
    --plotWidth 15 \
    --colors "#E41A1C" "#377EB8" "#4DAF4A" \
    --plotFileFormat png \
    --dpi 300

echo "=== TSS Metagene plots generated ==="
echo "  IP plot: $OUT_DIR/TSS_metagene_IP.png"
echo "  IN plot: $OUT_DIR/TSS_metagene_IN.png"

# ================================
# 8. Limpieza
# ================================
echo "Cleaning temporary files..."

rm -f "$VIS_DIR"/*.fixed.bedgraph
rm -f "$VIS_DIR"/*.bw
rm -f "$GTF_FIXED"
rm -f "$CHROM_SIZES"

echo "Cleanup done."
echo "=== Module completed successfully ==="