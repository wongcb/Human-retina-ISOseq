########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

cd /mnt/d/bioinformatic/isoseq/processed/from_Ray/proActiv

#Extract PBids for FL2 + CAGE + polyA motif
awk -F'\t' '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="isoform") iso_col=i
    if ($i=="within_cage_peak") cage_col=i
    if ($i=="polyA_motif") polya_col=i
  }

  if (!iso_col || !cage_col || !polya_col) {
    print "[ERROR] Cannot find required columns in classification file." > "/dev/stderr"
    print "isoform col=" iso_col ", within_cage_peak col=" cage_col ", polyA_motif col=" polya_col > "/dev/stderr"
    exit 1
  }
  next
}

($cage_col=="TRUE" && $polya_col!="" && $polya_col!="NA") {
  print $iso_col
}
' RetinaMerge_classification.filtered_lite_FL2_classification.txt \
  | sort -u > pbids_FL2_CAGE_polyA.txt

echo "Number of FL2+CAGE+polyA PBids:"
wc -l pbids_FL2_CAGE_polyA.txt
head pbids_FL2_CAGE_polyA.txt


#Extract molecule IDs from annotated info file
awk -F'\t' '
NR==FNR {
  keep[$1]=1
  next
}

FNR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="id") id_col=i
    if ($i=="pbid") pbid_col=i
  }

  if (!id_col || !pbid_col) {
    print "[ERROR] Cannot find required columns in annotated info file." > "/dev/stderr"
    print "id col=" id_col ", pbid col=" pbid_col > "/dev/stderr"
    exit 1
  }
  next
}

($pbid_col in keep) {
  print $id_col
}
' pbids_FL2_CAGE_polyA.txt RetinaMerge.annotated.info.csv \
  | sort -u > molecule_ids_FL2_CAGE_polyA.txt

echo "Number of FL2+CAGE+polyA molecule IDs:"
wc -l molecule_ids_FL2_CAGE_polyA.txt
head molecule_ids_FL2_CAGE_polyA.txt


#Split all bam files
cd RayMerged_bamfile

for bam in *.bam; do
  base=${bam%.bam}
  out="${base}_FL2_CAGE_polyA.bam"

  echo "Filtering $bam -> $out"
  samtools view -bh -N ../molecule_ids_FL2_CAGE_polyA.txt "$bam" > "$out"
  samtools index "$out"
done

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################