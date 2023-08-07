#!/usr/bin/env bash

input_files=${@}

echo -e "contig	dist.cM	weighted.LD	bin.count	bin.status"

for file in ${input_files}; do
  if [[ "${file}" =~ .*:.*$ ]]; then
    contig=$(basename ${file} | awk -F ":" '{print $2}')
  else
    contig="all"
  fi
  while read line; do
    dist=$(echo ${line} | awk '{print $1}')
    weighted_LD=$(echo ${line} | awk '{print $4}')
    bin_count=$(echo ${line} | awk '{print $5}')
    bin_status="NA"
    echo -e "${contig}\t${dist}\t${weighted_LD}\t${bin_count}\t${bin_status}"
  done < <(cat ${file} | grep -v "#" | sed 's/^[ \t]*//;s/[ \t]*$//' | tr -s " "  )
done