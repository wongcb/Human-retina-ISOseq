########Iso-seq profiling of novel isoform variants in human retina at single cell resolution############
########created by Luozixian Wang, Raymond Wong##########################################################
########CERA, UNIMELB, 20/06/2026########################################################################

cd /mnt/d/bioinformatic/isoseq/processed/from_Ray/ASprofile/FL2CAGEpolyAmotif

awk -F'\t' 'NR>1 && $41=="TRUE" && $44!="" && $44!="NA" && $44!="0" {print $1}' \
  RetinaMerge_classification.filtered_lite_classification.high_confident.txt \
  | sort -u > pbids_FL2_CAGE_polyA_strict.txt

awk '
NR==FNR {keep[$1]=1; next}
BEGIN {FS=OFS="\t"}
/^#/ {print; next}
{
  tid=""
  if (match($9, /transcript_id "([^"]+)"/)) {
    tid = substr($9, RSTART+15, RLENGTH-16)
  }
  if (tid in keep) print
}
' pbids_FL2_CAGE_polyA_strict.txt RetinaMerge.sorted.filtered_lite.gff \
> RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff

#check the output
echo "PBid count:"
wc -l pbids_FL2_CAGE_polyA_strict.txt

echo "Unique transcript count in subset gff:"
awk '
BEGIN {FS="\t"}
!/^#/ {
  if (match($9, /transcript_id "([^"]+)"/)) {
    print substr($9, RSTART+15, RLENGTH-16)
  }
}
' RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff | sort -u | wc -l

#ASprofile
/mnt/d/bioinformatic/isoseq/processed/from_Ray/ASprofile/ASprofile.b-1.0.4/extract-as RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff GRCh38.p14.genome.fa > AS_output_FL2CAGEpolyAmotif.txt

perl /mnt/d/bioinformatic/isoseq/processed/from_Ray/ASprofile/ASprofile.b-1.0.4/summarize_as.pl RetinaMerge.sorted.filtered_lite.FL2_CAGE_polyA_strict.gff AS_output_FL2CAGEpolyAmotif.txt -p Retina

####end of the session########################
####author: Luozixian Wang, Raymond Wong######
####20/06/2026################################