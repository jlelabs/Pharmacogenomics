#!/bin/bash

set -x

STARGAZER=${PWD}/Configuration/Stargazer_v1.0.8

while getopts ":i:v:o:b:r:" flag
do
    case "${flag}" in
        i) SAMPLE=${OPTARG};;
        v) VCF=${OPTARG};;
        b) BAM=${OPTARG};;
        o) PATHTOOUTPUT=${OPTARG};;
        r) REFERENCE=${OPTARG};;
        \?) valid=0
            echo "An invalid option has been entered: $OPTARG"
            exit 0
            ;;

        :)  valid=0
            echo "The additional argument for option $OPTARG was omitted."
            exit 0
            ;;

    esac
done

shift "$(( OPTIND - 1 ))"

if [ -z "$SAMPLE" ]; then
        echo 'Missing -i ID' >&2
        exit 1
fi

if [ -z "$VCF" ]; then
        echo 'Missing -v VCF file' >&2
        exit 1
fi

if [ -z "$BAM" ]; then
        echo 'Missing -b BAM file' >&2
        exit 1
fi

if [ -z "$PATHTOOUTPUT" ]; then
        echo 'Missing -o output directory' >&2
        exit 1
fi

if [ -z "$REFERENCE" ]; then
        echo 'Missing -r reference fasta' >&2
        exit 1
fi


echo $SAMPLE
echo $VCF
echo $BAM
echo $PATHTOOUTPUT
echo $REFERENCE


### folders
CONFIG=${PWD}/Configuration
OUTDIR=${PATHTOOUTPUT}/${SAMPLE}
PHARMACOGENOMICS=${OUTDIR}/${SAMPLE}_PHARMACOGENOMICS_REPORT
PHARMACOGENOMICS_ADDITIONAL=${OUTDIR}/${SAMPLE}_PHARMACOGENOMICS_SUPPLEMENTARY

mkdir -p ${OUTDIR}
mkdir -p ${PHARMACOGENOMICS}
mkdir -p ${PHARMACOGENOMICS_ADDITIONAL}

##### PHARMACOGENOMICS

{
  # FOR STARGAZER CONTROL RYR1
  gatk DepthOfCoverage \
    -R ${REFERENCE} \
    -I ${BAM} \
    -O ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_target_GDF_RYR1.gdf \
    -intervals ${CONFIG}/target_genes_RYR1.bed

  # convert GDF file from comma delimited to tab delimited
  tr ',' '\t' < ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_target_GDF_RYR1.gdf > \
    ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_target_GDF_RYR1_tab.gdf

  PHARMACOGENOMICS_TEMPDIR=${PHARMACOGENOMICS}/Temporary
  mkdir ${PHARMACOGENOMICS_TEMPDIR}

  # structural variant with gdf for cyp2d6
  #ryr1
  python ${STARGAZER}/stargazer.py genotype \
    -o ${PHARMACOGENOMICS_TEMPDIR}/${SAMPLE}.cyp2d6 \
    -d wgs \
    -t cyp2d6 \
    -c ryr1 \
    --vcf ${VCF} \
    --gdf ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_target_GDF_RYR1_tab.gdf

  for gene in cyp2c9 cyp2c19 cyp4f2 vkorc1 dpyd g6pd ifnl3 nudt15 slco1b1 tpmt ugt1a1 cyp3a5 cacna1s ryr1 cyp2b6 cftr; do
    python ${STARGAZER}/stargazer.py genotype \
      -o ${PHARMACOGENOMICS_TEMPDIR}/${SAMPLE}.${gene} \
      -d wgs \
      -t ${gene} \
      --vcf ${VCF}
  done

  echo "" > ${PHARMACOGENOMICS_TEMPDIR}/${SAMPLE}_temp.txt
  header=""

  for g in ${PHARMACOGENOMICS_TEMPDIR}/*.stargazer-genotype.txt; do
    file_name=$(basename $g)
    gene=$(cut -d'.' -f2 <<<"$file_name")
    awk -v var="$gene" 'NR==1 {printf("%s\t%s\n", $0, "gene")}  NR>1 {printf("%s\t%s\n", $0, var) }' $g > ${PHARMACOGENOMICS_TEMPDIR}/${gene}.new.txt
    awk 'NR==2' ${PHARMACOGENOMICS_TEMPDIR}/${gene}.new.txt >> ${PHARMACOGENOMICS_TEMPDIR}/${SAMPLE}_temp.txt
    header=`awk 'NR==1' ${PHARMACOGENOMICS_TEMPDIR}/${gene}.new.txt`
  done

  echo "${header}" > ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_final_stargazer_output.txt
  cat ${PHARMACOGENOMICS_TEMPDIR}/${SAMPLE}_temp.txt >> ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_final_stargazer_output.txt

  # checking for rs12777823 for warfarin
  bcftools view --targets "chr10:96405502" -H ${VCF} > ${PHARMACOGENOMICS_ADDITIONAL}/${SAMPLE}_rsID_Warfarin.txt

  python calling_pharma.py ${SAMPLE} ${PHARMACOGENOMICS} ${PHARMACOGENOMICS_ADDITIONAL} ${HLA}
  rm -r ${PHARMACOGENOMICS_TEMPDIR}

} 2>&1 | tee ${PHARMACOGENOMICS}/Pharmacogenomics.log
