#!/bin/bash

# ========== path==========
HUMAN_ID_FILE="/home/bsb/mywork/data/middlefile/human_uniprot_ids.txt"
MOUSE_ID_FILE="/home/bsb/mywork/data/middlefile/mouse_uniprot_ids.txt"
HUMAN_SRC_DIR="/home/bsb/mywork/rawdata/structure/human"
MOUSE_SRC_DIR="/home/bsb/mywork/rawdata/structure/mouse"
HUMAN_DST_DIR="/home/bsb/mywork/data/middlefile/human"
MOUSE_DST_DIR="/home/bsb/mywork/data/middlefile/mouse"
# 
echo "Copying human structures..."
while IFS= read -r id; do
  if [ -n "$id" ]; then
    # match: AF-{id}-F1-model_v6.pdb 
    find "$HUMAN_SRC_DIR" -maxdepth 1 -name "AF-${id}-F1-model_v6.pdb" -exec cp {} "$HUMAN_DST_DIR/" \; 2>/dev/null
  fi
done < "$HUMAN_ID_FILE"

# 
echo "Copying mouse structures..."
while IFS= read -r id; do
  if [ -n "$id" ]; then
    find "$MOUSE_SRC_DIR" -maxdepth 1 -name "AF-${id}-F1-model_v6.pdb" -exec cp {} "$MOUSE_DST_DIR/" \; 2>/dev/null
  fi
done < "$MOUSE_ID_FILE"

# ========== sum ==========
human_count=$(ls "$HUMAN_DST_DIR" 2>/dev/null | wc -l)
mouse_count=$(ls "$MOUSE_DST_DIR" 2>/dev/null | wc -l)

echo "✅ Done!"
echo "Human structures copied: $human_count"
echo "Mouse structures copied: $mouse_count"
cp /home/bsb/mywork/rawdata/structure/missing_structure/* /home/bsb/mywork/data/middlefile/mouse
