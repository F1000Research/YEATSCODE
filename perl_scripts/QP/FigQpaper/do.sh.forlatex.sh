#!/bin/tcsh


latex document.tex
bibtex document
latex document.tex
latex document.tex

#dvips -Ppdf -o document.ps document.dvi
#ps2pdf document.ps document.pdf


