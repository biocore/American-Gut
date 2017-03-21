#!/bin/bash

# Grab any sample accessions that do not exist

study_accession=ERP012803
count=-1

# Grab all the secondary sample accession numbers and file links from EBI
curl -s "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=${study_accession}&result=read_run&fields=secondary_sample_accession,submitted_ftp" | grep -v "^secondary_sample_accession" > ${study_accession}.details.txt

# Each line of the file `${study_accession}.details.txt` contains the secondary
# sample accession followed by the link to download the file
# This `awk` command splits such line in 2
for fq in `awk '{print $1, $2}' ${study_accession}.details.txt`
do
    ((count++))
    if [[ $(( count % 2)) -eq 0 ]]
    then
        # Even case, we have the secondary sample accession number
        id=$fq
        current_path=${study_accession}/${id}
        current_base=${current_path}/${id}

        # Check if we already have the metadata of the currrent sample
        if [ -d "${current_path}" ]; then
            continue
        fi
        echo "Fetching ${id}..."

        # Retrieving the sample metadata
        mkdir -p ${current_path}
        curl -s "http://www.ebi.ac.uk/ena/data/view/${id}&display=xml" > ${current_base}.txt &
    else
        # Odd case, we have the link to the sequence file
        # Check if we already have the sequence file for the current sample,
        # note that current_base gets defined in the previous iteration of the
        # for loop
        if [ -e "${current_base}.fna" ]; then
            continue
        fi

        # Retrieve the FASTQ sequence file and transform it to FASTA
        # sed from http://stackoverflow.com/a/10359425/19741
        curl -s $fq | zcat | sed -n '1~4s/^@/>/p;2~4p' > ${current_base}.fna &
    fi

    # Only doing 10 downloads (5 samples) in parallel, so wait until the
    # previous processes are done before continuing downloading data
    if [[ $((count % 10)) -eq 0 ]]
    then
        wait
    fi
done

# Wait until all downloads are completed
wait
