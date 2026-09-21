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

.PHONY: all figures refs paper paper-draft test audit check-artefacts clean distclean symbolic

all: figures paper

# The three model modules keep their old directory names; the paper no longer
# uses the letters. model_C is the dynamic model and its equilibrium
# reduction, model_A is the local identification model, model_B is the
# alternative hypothesis with a ceiling (one remark in the paper now).
figures:
	@echo ">> the dynamic model and its equilibrium reduction"; $(PY) src/model_C/run.py
	@echo ">> the local identification model";                  $(PY) src/model_A/run.py
	@echo ">> the ceiling hypothesis";                          $(PY) src/model_B/run.py
	@echo ">> the design (slow: the decision characteristic is simulated)"; $(PY) src/design/run.py
	@echo ">> the scope tables, whose convention lives in the audit script"; $(PY) audit/n4_theta_scope.py
	@echo ">> the decision table (slow: 12 seeds x 400 replicates)";        $(PY) audit/d13_decision_therapy.py
	@echo ">> the admissibility checks and the dose sensitivity";           $(PY) audit/h1_admissibility.py
	@echo ">> the dose calibration table (slow)";                           $(PY) audit/k4_admissibility.py
	@echo ">> the K_m sweep and the profile interval, which the text quotes"
	$(PY) audit/n1_km_cost.py
	$(PY) audit/n6_profile_vs_crb.py
	@echo ">> the transient's worth, which Sections 4.1 and 6.2 quote"
	$(PY) audit/t1_transient_worth.py
	@echo ">> the band against the VALUE of k_e, which Section 6.2 quotes"
	$(PY) audit/t5_band_vs_ke.py
	@echo ">> the cohort coverage with the increment CV drawn (Section 5.5)"
	$(PY) audit/t6_cohort_with_kappa.py
	@echo ">> what the structure test costs (Section 5.1)"
	$(PY) audit/t7_structure_transient.py
	@echo ">> the phase lattice and the lap condition (Remark 4)"
	$(PY) audit/t11_phase_lattice.py
	@echo ">> the forward therapy simulation, whose run (b) Section 4.2 quotes"
	$(PY) audit/t13_therapy_forward.py
	@echo ">> the design figures: phase clock, visits, calendar, conditioning"
	$(PY) audit/g1_phase_clock.py
	$(PY) audit/g2_visits_vs_sd.py
	$(PY) audit/g3_calendar.py
	$(PY) audit/g5_conditions.py

refs:
	$(PY) tools/fetch_refs.py
	$(PY) tools/claims_to_tex.py

# The paper and the supplement are separate documents. The supplement imports
# the paper's labels with xr, so the paper MUST be built first; the paper cites
# the supplement by section number in prose, so there is no cycle to iterate.
paper: refs
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && biber main >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode supplement.tex >/dev/null
	cd paper && biber supplement >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode supplement.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode supplement.tex >/dev/null
	@echo ">> built $(PAPER).pdf and paper/supplement.pdf"

symbolic:
	@echo "Structural identifiability requires Julia:"
	@echo "  julia -e 'using Pkg; Pkg.add(\"StructuralIdentifiability\")'"
	@echo "  julia symbolic/julia/run_identifiability.jl"

test:
	$(PY) -m pytest tests/ -q

# Every committed table, regenerated and compared against what is on disk.
# `make test` runs the same check; this target is for running it alone.
# --working-tree instead compares all of results/ against HEAD, figures
# included, and is what CI runs after `make figures`.
check-artefacts:
	$(PY) tools/check_artefacts.py

# Non-zero while any claim is contradicted or unsourced, or any identifier
# fails to resolve. Today it fails on 14 unresolved claims, which is the
# correct state of the audit and not a broken build.
# The same document with the claim-register appendix appended. The register
# is a working artefact, so it is off in the default build.
paper-draft: refs
	cd paper && $(LATEX) -interaction=nonstopmode -jobname=main \
	  "\\PassOptionsToClass{draft}{article}\\input{main}" >/dev/null
	cd paper && biber main >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode -jobname=main \
	  "\\PassOptionsToClass{draft}{article}\\input{main}" >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode -jobname=main \
	  "\\PassOptionsToClass{draft}{article}\\input{main}" >/dev/null
	@echo ">> built $(PAPER).pdf with the claim register"

audit:
	$(PY) tools/fetch_refs.py --check
	$(PY) tools/claims_to_tex.py --check

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc paper/*.fls paper/*.fdb_latexmk
	rm -f paper/*.bbl paper/*.bcf paper/*.blg paper/*.run.xml
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

distclean: clean
	rm -f results/figures/* results/tables/* paper/main.pdf
