#Plotting allele frequency of LP-associated rs4988235 with % of vegetarians in Indian populations from 1000 Genomes 30x GRCh38 collection using 1991 data
import numpy as np
import matplotlib.pyplot as plt
#Data is from Nutrition journal and vegetarian corresponds to lacto-vegetarian in Table 2
af = [0.140777,0.236301,0.035088,0.056075,0.057252]
veg = [55.6,52.3,4.5,4.7,1.4]
ethnicities = ['Gujarati','Punjabi','Tamil','Telugu','Bengali']

correlation_matrix = np.corrcoef(af, veg)
r_squared = correlation_matrix[0, 1] ** 2

plt.scatter(af, veg, color='blue', edgecolors='black')

for i, ethnicity in enumerate(ethnicities):
    plt.annotate(
        ethnicity,                      # The text label
        (af[i], veg[i]), # The point (x, y) to label
        textcoords="offset points", # How to position the text
        xytext=(0, 10),             # Distance from text to point (x, y)
        ha='center'                 # Horizontal alignment (center, left, or right)
    )

plt.title("LP allele frequency against Vegetarian proportion")
plt.xlabel("rs4988235 allele frequency")
plt.ylabel("Percentage of vegetarians")

plt.show()

print(r_squared)





#Download the full GRCh37 reference
cd ~/thesis_data/required_files
wget -b -c -O human_g1k_v37.fasta.gz \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

#Integrity check (important — file misreads this multi-stream gzip)
gzip -t human_g1k_v37.fasta.gz; echo "Exit: $?"
# Exit 2 with "trailing garbage ignored" is acceptable; decompresses fine.

#Decompress and extract chromosome 2
gunzip human_g1k_v37.fasta.gz
# If gunzip refuses due to warning:
zcat human_g1k_v37.fasta.gz > human_g1k_v37.fasta 2>/dev/null

# samtools failed on this system due to libncurses.so.5 — awk fallback:
awk '/^>2 /{flag=1; print; next} /^>/{flag=0} flag' \
    human_g1k_v37.fasta > human_g1k_v37_chr2.fa

#Rename header to match downstream tools
sed -i '1s/^>2.*$/>2/' human_g1k_v37_chr2.fa
head -1 human_g1k_v37_chr2.fa   # Expected: >2

#Verify
grep -v "^>" human_g1k_v37_chr2.fa | tr -d '\n' | wc -c   # 243199373
seqkit stats human_g1k_v37_chr2.fa

#Download the 1000 Genomes ancestral sequence
cd ~/thesis_data/required_files
wget -O human_ancestor_2.fa \
  "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_2.fa"

# Alternative (EBI FTP):
wget -O human_ancestor_2.fa \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_2.fa"

#Diagnose the "dot convention"
grep -v "^>" human_ancestor_2.fa | tr -d '\n' | \
  awk '{n_dot=gsub(/\./,"."); n_other=length($0)-n_dot;
        print "Dots: " n_dot, "  Non-dot: " n_other}'

grep -v "^>" human_ancestor_2.fa | tr -d '\n' | \
  grep -o "[ACGTNacgtn]" | sort | uniq -c

#Reconstruct the ancestral FASTA (fill dots from reference, uppercase everything)
cd ~/thesis_data/required_files

python3 - <<'EOF'
def read_fasta(path):
    with open(path) as f:
        header = f.readline().strip()
        seq = ''.join(line.strip() for line in f)
    return header, seq

_, ref = read_fasta("human_g1k_v37_chr2.fa")
_, anc = read_fasta("human_ancestor_2.fa")

print(f"Reference length: {len(ref)}")
print(f"Ancestral length: {len(anc)}")

if len(ref) != len(anc):
    raise SystemExit("Length mismatch — check both files cover chromosome 2 (GRCh37).")

merged = []
for a, r in zip(anc, ref):
    if a == ".":
        merged.append(r.upper())     # dot means "matches reference"
    elif a in "ACGTN":
        merged.append(a)             # explicit ancestral base
    elif a in "acgtn":
        merged.append(a.upper())     # soft-masked base -> uppercase for Relate
    else:
        merged.append("N")           # anything else -> unknown

merged = ''.join(merged)


with open("human_ancestor_2_reconstructed.fa", "w") as out:
    out.write(">2\n")
    for i in range(0, len(merged), 60):
        out.write(merged[i:i+60] + "\n")

print(f"Wrote human_ancestor_2_reconstructed.fa (length {len(merged)})")
EOF

#Rename header and verify
sed -i '1s/^>.*$/>2/' human_ancestor_2_reconstructed.fa

head -1 human_ancestor_2_reconstructed.fa                    # >2
grep -v "^>" human_ancestor_2_reconstructed.fa | tr -d '\n' | wc -c   # 243199373

# Confirm no dots and no lowercase remain
grep -v "^>" human_ancestor_2_reconstructed.fa | tr -d '\n' | \
  awk '{print "Dots: " gsub(/\./,"."), " Lowercase: " gsub(/[acgtn]/,"a")}'

#Focal-site validation (before Phase 1)
# Base at rs4988235 (must be G — the ancestral allele, per VCF INFO AA=G)
grep -v "^>" human_ancestor_2_reconstructed.fa | tr -d '\n' | cut -c136608646

# N-content in a 20 kb window around the SNP
grep -v "^>" human_ancestor_2_reconstructed.fa | tr -d '\n' | cut -c136600000-136620000 | \
  awk '{n=gsub(/N/,"N"); print "N count in 20 kb: " n " (" n/20000*100 "%)"}'
# Observed: 28 N (0.14%) — excellent

#Download the Phase 1 strict mask (FASTA, per-chromosome)
cd ~/thesis_data/required_files
wget -O chr2.strictMask.fasta.gz \
  "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/supporting/accessible_genome_masks/chr2.strictMask.fasta.gz"
gunzip chr2.strictMask.fasta.gz

#Or download the Phase 3 strict mask (BED, whole-genome)
wget -O 20141020.strict_mask.whole_genome.bed \
  "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed"

#Convert the six-code Phase 1 mask to binary P/N
cd ~/thesis_data/required_files

awk 'BEGIN{FS=""}
     /^>/ {print ">2"; next}
     {for(i=1;i<=NF;i++){
         c=$i;
         if(c=="P") printf "P";
         else printf "N";
     }
     printf "\n"}' chr2.strictMask.fasta > 20140520.chr2.strict_mask.fasta

#Verify
head -1 20140520.chr2.strict_mask.fasta                          # >2
grep -v "^>" 20140520.chr2.strict_mask.fasta | tr -d '\n' | wc -c   # 243199373

# P count should match the header annotation (~179021843)
grep -v "^>" 20140520.chr2.strict_mask.fasta | tr -d '\n' | \
  awk '{n=gsub(/P/,"P"); print "P count: " n}'

# Focal-site accessibility
grep -v "^>" 20140520.chr2.strict_mask.fasta | tr -d '\n' | cut -c136608646
# Expected: P

#Note on inspecting a FASTA (line-wrapping pitfall)
grep -v "^>" file.fa | head -5 | tr -d '\n' | cut -c1-300

#Download (OMNI-interpolated, GRCh37)
cd ~/thesis_data/required_files
wget -O genetic_map_chr2.txt \
  "https://raw.githubusercontent.com/DrKLWilliams/1000-genomes-genetic-maps/master/interpolated_OMNI/chr2.txt"

#Verify format
head -5 genetic_map_chr2.txt
# Observed (correct):
#   position    rate    GeneticMap
#   12994       0       0
#   15491       0.339607529034842  0.000848

sed -i '1s/.*/position\tCOMBINED_rate(cM\/Mb)\tGenetic_Map(cM)/' genetic_map_chr2.txt
head -1 genetic_map_chr2.txt

#Verify coverage of the focal site
awk '$1 <= 136608646' genetic_map_chr2.txt | tail -1
awk '$1 >= 136608646' genetic_map_chr2.txt | head -1
tail -1 genetic_map_chr2.txt


#Download 1000 Genomes Phase 3 chromosome 2
mkdir -p ~/thesis_data/raw_vcf
cd ~/thesis_data/raw_vcf

wget -O ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

wget -O ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"

#Download the sample panel
wget -O integrated_call_samples_v3.20130502.ALL.panel \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

#Create per-population sample lists
for pop in BEB GIH ITU PJL STU; do
    grep "$pop" integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > ${pop}_samples.txt
done

#Subset VCF to SAS samples and the LCT/MCM6 region
grep -E "BEB|GIH|ITU|PJL|STU" integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > SAS_samples.txt

bcftools view -S SAS_samples.txt \
  -r 2:136000000-137000000 \
  -Oz -o chr2_SAS_LCT.vcf.gz \
  ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools index -t chr2_SAS_LCT.vcf.gz

#Split into per-population VCFs with Companion-compatible names
for pop in BEB GIH ITU PJL STU; do
    bcftools view -S ${pop}_samples.txt \
      -Oz -o chr2_${pop}_LCT.vcf.gz \
      chr2_SAS_LCT.vcf.gz
    bcftools index -t chr2_${pop}_LCT.vcf.gz
done

#Annotate rs4988235 (the 1000G VCF omits rsIDs)
cd ~/thesis_data/raw_vcf

# Tab-separated, with header
printf "#CHROM\tPOS\tID\n2\t136608646\trs4988235\n" > rs4988235_annot.txt

# Compress and index the annotation file (bcftools annotate requires this)
bgzip rs4988235_annot.txt
tabix -s1 -b2 -e2 rs4988235_annot.txt.gz

for pop in BEB GIH ITU PJL STU; do
    bcftools annotate \
      -a rs4988235_annot.txt.gz \
      -c CHROM,POS,ID \
      -Oz -o ${pop}_LCT_annotated_chr2.vcf.gz \
      chr2_${pop}_LCT.vcf.gz
    bcftools index -t ${pop}_LCT_annotated_chr2.vcf.gz
done

#Create per-population .poplabels files
cd ~/thesis_data/poplabels
for pop in BEB GIH ITU PJL STU; do
    echo -e "sample\tpopulation\tgroup\tsex" > ${pop}.poplabels
    grep "$pop" ~/thesis_data/raw_vcf/integrated_call_samples_v3.20130502.ALL.panel | \
      awk -v p=$pop '{print $1"\t"p"\tSAS\tNA"}' >> ${pop}.poplabels
done

#6. Pre-Phase-1 Sanity Checks (Reproducible)
#6.1 VCF sanity
for pop in BEB GIH ITU PJL STU; do
    f=~/thesis_data/raw_vcf/${pop}_LCT_annotated_chr2.vcf.gz
    echo "=== $pop ==="
    bcftools view -h "$f" | grep "^##contig=<ID=2,"           # chromosome 2 defined
    bcftools view "$f" | grep -v "^#" | cut -f1 | sort -u     # variants on chr 2 only
    bcftools query -l "$f" | wc -l                            # sample count
    bcftools view -r 2:136608646 "$f" | grep -v "^#" | cut -f1-5
done

#6.2 Derived allele frequency per population
for pop in BEB GIH ITU PJL STU; do
    f=~/thesis_data/raw_vcf/${pop}_LCT_annotated_chr2.vcf.gz
    echo -n "$pop  DAF = "
    bcftools view -H -r 2:136608646 "$f" | \
      awk '{for(i=10;i<=NF;i++){split($i,a,"|"); if(a[1]==1) d++; if(a[2]==1) d++; n+=2}}
           END{printf "%.4f  (%d/%d)\n", d/n, d, n}'
done

#6.3 Focal-site cross-consistency
for pop in BEB GIH ITU PJL STU; do
    f=~/thesis_data/raw_vcf/${pop}_LCT_annotated_chr2.vcf.gz
    echo -n "$pop REF: "
    bcftools view -r 2:136608646 "$f" | grep -v "^#" | cut -f4
done
echo -n "Ancestral: "
grep -v "^>" ~/companion_clues2/CLUES2-Companion/required_files/ancestor/homo_sapiens_ancestor_chr2.fa | \
    tr -d '\n' | cut -c136608646

#6.4 Mask accessibility at focal site
grep -v "^>" ~/companion_clues2/CLUES2-Companion/required_files/mask/PilotMask_chr2.fasta | \
    tr -d '\n' | cut -c136608646      # Expected: P

#6.5 Recombination map covers focal site
awk '$1 <= 136608646' ~/companion_clues2/CLUES2-Companion/required_files/map/genetic_map_chr2.txt | tail -1
awk '$1 >= 136608646' ~/companion_clues2/CLUES2-Companion/required_files/map/genetic_map_chr2.txt | head -1

#Install the Auxiliary Files Where the Companion Expects Them
cd ~/companion_clues2/CLUES2-Companion

mkdir -p required_files/ancestor
mkdir -p required_files/mask
mkdir -p required_files/map

cp ~/thesis_data/required_files/human_ancestor_2_reconstructed.fa \
   required_files/ancestor/homo_sapiens_ancestor_chr2.fa

cp ~/thesis_data/required_files/20140520.chr2.strict_mask.fasta \
   required_files/mask/PilotMask_chr2.fasta

cp ~/thesis_data/required_files/genetic_map_chr2.txt \
   required_files/map/genetic_map_chr2.txt

find ~/companion_clues2/CLUES2-Companion/required_files -type f | sort

#Launching Phase 1
cd ~/companion_clues2/CLUES2-Companion
./CLUES2Companion.sh

#Verify Phase 1 output
ls -lh ~/companion_clues2/CLUES2-Companion/output_CLUES2Companion-Relate/phase1/PJL_MCM6_chr2/

grep rs4988235 ~/companion_clues2/CLUES2-Companion/output_CLUES2Companion-Relate/phase1/PJL_MCM6_chr2/*SNPs*
grep rs4988235 ~/companion_clues2/CLUES2-Companion/output_CLUES2Companion-Relate/phase1/PJL_MCM6_chr2/*Frequency*
