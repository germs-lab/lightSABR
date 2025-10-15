#!/usr/bin/env bash
set -Eeuo pipefail


# Usage: 
# ./blast_nt_remote.sh ASV.fasta blast_remote.tsv


# Remote BLAST against nt with an Entrez query excluding uninformative titles.

FASTA="$1"
OUT="$2"
THREADS="${3:-2}"

if [[ -z "$FASTA" || -z "$OUT" ]]; then
  echo "Usage: $0 asv.fasta blast_remote.tsv [threads]" >&2
  exit 1
fi

ENTREZ_QUERY='"Bacteria"[Organism] NOT (uncultured[Title] OR environmental[Title] OR metagenome[Title] OR clone[Title] OR unidentified[Title] OR "sp." [Title])'

blastn \
  -task megablast \
  -db nt \
  -remote \
  -query "$FASTA" \
  -evalue 1e-20 \
  -perc_identity 99 \
  -qcov_hsp_perc 90 \
  -max_target_seqs 25 \
  -max_hsps 1 \
  -dust no \
  -entrez_query "$ENTREZ_QUERY" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids ssciname stitle' \
  -out "$OUT"

# BLAST 2.13.0+
# USAGE
#   blastn [-h] [-help] [-import_search_strategy filename]
#     [-export_search_strategy filename] [-task task_name] [-db database_name]
#     [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
#     [-negative_gilist filename] [-negative_seqidlist filename]
#     [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
#     [-negative_taxidlist filename] [-entrez_query entrez_query]
#     [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
#     [-subject subject_input_file] [-subject_loc range] [-query input_file]
#     [-out output_file] [-evalue evalue] [-word_size int_value]
#     [-gapopen open_penalty] [-gapextend extend_penalty]
#     [-perc_identity float_value] [-qcov_hsp_perc float_value]
#     [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
#     [-xdrop_gap_final float_value] [-searchsp int_value]
#     [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
#     [-min_raw_gapped_score int_value] [-template_type type]
#     [-template_length int_value] [-dust DUST_options]
#     [-filtering_db filtering_database]
#     [-window_masker_taxid window_masker_taxid]
#     [-window_masker_db window_masker_db] [-soft_masking soft_masking]
#     [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
#     [-best_hit_score_edge float_value] [-subject_besthit]
#     [-window_size int_value] [-off_diagonal_range int_value]
#     [-use_index boolean] [-index_name string] [-lcase_masking]
#     [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
#     [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
#     [-line_length line_length] [-html] [-sorthits sort_hits]
#     [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
#     [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

