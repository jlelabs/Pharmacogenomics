#!/bin/bash

set -x

while getopts ":r:" flag
do
    case "${flag}" in
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

if [ -z "$REFERENCE" ]; then
        echo 'Missing -r reference fasta' >&2
        exit 1
fi

REFERENCE_DIR=$(dirname ${REFERENCE})
REFERENCE_NAME=$(basename ${REFERENCE})
REFERENCE_BASE="${REFERENCE_NAME%.*}"

samtools faidx ${REFERENCE}
gatk CreateSequenceDictionary -R ${REFERENCE} -O ${REFERENCE_DIR}/${REFERENCE_BASE}.dict
