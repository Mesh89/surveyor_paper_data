VCF=$1
LR_FQ=$2
OUTDIR=$3
OUT_PREFIX=$4
THREADS=$5

mkdir -p $OUTDIR
if ~/arm/bin/bcftools view -h "$VCF" | grep -q '^##FORMAT=<ID=FT,'; then
    # FT defined: require FT to be PASS
    echo "FT defined: require FT to be PASS"
    FILTER_EXPR="SVTYPE=='INV' && FT=='PASS'"
else
    # FT not defined: ignore the FT condition
    echo "FT not defined: ignore the FT condition"
    FILTER_EXPR="SVTYPE=='INV'"
fi
 bcftools view -i "$FILTER_EXPR" --min-ac=1 -f PASS,. "$VCF" -Oz -o "$OUTDIR/$OUT_PREFIX.vcf.gz"
./validate.sh $OUTDIR/$OUT_PREFIX.vcf.gz $LR_FQ $OUTDIR/ $THREADS
