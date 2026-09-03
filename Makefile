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

.PHONY: all figures paper test clean distclean symbolic

all: figures paper

figures:
	@echo ">> Model A"; $(PY) src/model_A/run.py
	@echo ">> Model B"; $(PY) src/model_B/run.py
	@echo ">> Model C"; $(PY) src/model_C/run.py

paper:
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo ">> built $(PAPER).pdf"

symbolic:
	@echo "Structural identifiability requires Julia:"
	@echo "  julia -e 'using Pkg; Pkg.add(\"StructuralIdentifiability\")'"
	@echo "  julia symbolic/julia/run_identifiability.jl"

test:
	$(PY) -m pytest tests/ -q

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

distclean: clean
	rm -f results/figures/* results/tables/* paper/main.pdf
