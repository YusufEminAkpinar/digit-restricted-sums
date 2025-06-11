FILE = rapport3
TC = pdflatex
BIB = bibtex

all: compile

compile: $(FILE).tex
	$(TC) -shell-escape $(FILE).tex

bib: $(FILE).tex
	$(TC) -shell-escape $(FILE).tex
	$(BIB) $(FILE) || true # Ensure that make command continue
	$(TC) -shell-escape $(FILE).tex
	$(TC) -shell-escape $(FILE).tex


clean:
	rm *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.xml *.bcf *_latexmk *.fls 
