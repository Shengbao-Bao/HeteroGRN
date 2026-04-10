#!/bin/bash

HUMAN_DIR="/home/bsb/mywork/data/middlefile/human"
MOUSE_DIR="/home/bsb/mywork/data/middlefile/mouse"
HUMAN_RESULT="/home/bsb/mywork/data/middlefile/humanresult"
MOUSE_RESULT="/home/bsb/mywork/data/middlefile/mouseresult"

mkdir -p "$HUMAN_RESULT" "$MOUSE_RESULT"

run_tf_similarity() {
  local input_dir="$1"
  local output_file="$2"
  local tmp_dir="$3"
  local species_name="$4"

  echo "begin $species_name ..."

  mkdir -p "$tmp_dir"

  
  foldseek easy-search \
    --gpu 1 \
    --threads $(nproc) \
    --exhaustive-search 1 \
    --chain-name-mode 0 \
    --prefilter-mode 2 \
    --tmscore-threshold 0 \
    --lddt-threshold 0 \
    --min-aln-len 1 \
    -e 1000 \
    --format-output "query,target,alntmscore,prob,lddt,qtmscore,ttmscore,evalue" \
    "$input_dir" "$input_dir" "$output_file" "$tmp_dir"

  rm -rf "$tmp_dir"
  echo "✅ $species_name finished！"
}

run_tf_similarity "$HUMAN_DIR" "$HUMAN_RESULT/tf_similarity.tsv" "$HUMAN_RESULT/tmp" "human"
run_tf_similarity "$MOUSE_DIR" "$MOUSE_RESULT/tf_similarity.tsv" "$MOUSE_RESULT/tmp" "mouse"

echo "🎉 all finished！"
