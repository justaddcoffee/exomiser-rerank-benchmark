#!/usr/bin/env bash
# Fetch the Phenopacket Store bundle into data/ppkts/ (needs gh + curl + unzip).
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
