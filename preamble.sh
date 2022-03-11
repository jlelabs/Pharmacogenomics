#!/bin/bash

set -x

samtools faidx ${REFERENCE}
gatk CreateSequenceDictionary -R ${REFERENCE} -O ${REFERENCE_DIR}/${REFERENCE_BASE}.dict
