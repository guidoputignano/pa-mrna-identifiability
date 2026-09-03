# ===================================================================
#  propionate-identifiability
#  make figures   regenerate every figure and table from src/
#  make paper     build paper/main.pdf
#  make all       both
#  make test      run the consistency checks
#  make clean     remove build artefacts (figures are kept)
# ===================================================================

PY      ?= python
LATEX   ?= pdflatex
PAPER    = paper/main

.PHONY: all figures refs paper test audit clean distclean symbolic

all: figures paper

figures:
	@echo ">> Model A"; $(PY) src/model_A/run.py
	@echo ">> Model B"; $(PY) src/model_B/run.py
	@echo ">> Model C"; $(PY) src/model_C/run.py

refs:
	$(PY) tools/fetch_refs.py
	$(PY) tools/claims_to_tex.py

paper: refs
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && biber main >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo ">> built $(PAPER).pdf"

symbolic:
	@echo "Structural identifiability requires Julia:"
	@echo "  julia -e 'using Pkg; Pkg.add(\"StructuralIdentifiability\")'"
	@echo "  julia symbolic/julia/run_identifiability.jl"

test:
	$(PY) -m pytest tests/ -q

# Non-zero while any claim is contradicted or unsourced, or any identifier
# fails to resolve. Today it fails on 14 unresolved claims, which is the
# correct state of the audit and not a broken build.
audit:
	$(PY) tools/fetch_refs.py --check
	$(PY) tools/claims_to_tex.py --check

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk
	rm -f paper/*.bbl paper/*.bcf paper/*.blg paper/*.run.xml
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

distclean: clean
	rm -f results/figures/* results/tables/* paper/main.pdf
