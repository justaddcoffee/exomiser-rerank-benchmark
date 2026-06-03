#!/usr/bin/env bash
# Fetch the Phenopacket Store bundle into data/ppkts/, plus the HPO annotation files
# used by benchmark.synthesize (needs gh + curl + unzip).
set -euo pipefail
REPO=monarch-initiative/phenopacket-store
TAG=$(gh api "repos/$REPO/releases/latest" --jq .tag_name)
URL=$(gh api "repos/$REPO/releases/latest" \
  --jq '.assets[] | select(.name=="all_phenopackets.zip") | .browser_download_url')
mkdir -p data
curl -sL -o data/all_phenopackets.zip "$URL"
rm -rf data/ppkts && mkdir -p data/ppkts
unzip -q data/all_phenopackets.zip -d data/ppkts
echo "Phenopacket Store $TAG -> data/ppkts/ ($(find data/ppkts -name '*.json' | wc -l | tr -d ' ') phenopackets)"

# HPO annotation files for synthetic-case generation (benchmark.synthesize).
echo "Fetching HPOA + genes_to_disease + hp.obo for synthetic generation..."
curl -sL -o data/phenotype.hpoa       http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
curl -sL -o data/genes_to_disease.txt http://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt
curl -sL -o data/hp.obo               http://purl.obolibrary.org/obo/hp.obo
echo "  $(grep -m1 '#version:' data/phenotype.hpoa)"
