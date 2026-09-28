#  The change report: numbers, figures, then the document.
#    make            everything
#    make numbers    numbers.tex and provenance.csv only (seconds)
#    make figures    the four figures (about a minute: two worktree builds)
#    make section S=03_genetic   compile one section on its own
PY ?= python
TEX = pdflatex -interaction=nonstopmode -halt-on-error

all: report.pdf

numbers: numbers.tex
numbers.tex provenance.csv: make_numbers.py model_abc.py
	$(PY) make_numbers.py

figures: fig/phaseclock_0920.pdf fig/fig5a_abc.pdf
fig/phaseclock_0920.pdf fig/phaseclock_0927.pdf: fig_phaseclock.py
	$(PY) fig_phaseclock.py
fig/fig5a_abc.pdf fig/fig5b_fibre.pdf: fig_abc.py model_abc.py
	$(PY) fig_abc.py

report.pdf: report.tex sections/*.tex numbers.tex figures
	$(TEX) report.tex >/dev/null && $(TEX) report.tex | grep -E "Output written|Warning.*undefined" || true
	@grep -c "rep@num" numbers.tex >/dev/null

#  one section alone, so a failure in it costs nothing elsewhere
section: numbers.tex
	$(PY) one_section.py $(S)
	$(TEX) -jobname=_one_$(S) _one_$(S).tex >/dev/null \
	  && $(PY) -c "import pymupdf;print('sections/$(S).tex compiles,', pymupdf.open('_one_$(S).pdf').page_count, 'page(s)')" \
	  || grep -E "^!|^l\.[0-9]" _one_$(S).log
	@rm -f _one_$(S).*

clean:
	rm -f *.aux *.log *.out _one_*

.PHONY: all numbers figures section clean
