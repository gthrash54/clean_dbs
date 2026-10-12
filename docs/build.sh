#!/usr/bin/env bash
#builds every PDF in docs/ from its .tex (needs a LaTeX install with latexmk; same fonts as the CV)
set -e
cd "$(dirname "$0")"
for tex in user_guide.tex cheaha.tex; do
    latexmk -pdf -interaction=nonstopmode -halt-on-error -quiet "$tex"
done
latexmk -c -quiet user_guide.tex cheaha.tex #clean up the helper files, keep the PDFs
rm -f *.listing
