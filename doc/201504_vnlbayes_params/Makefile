TARGET = slides_vnlb_prms

finalpdf: $(TARGET).tex $(TARGET).pdf

zip: $(TARGET).zip

$(TARGET).zip: $(TARGET).pdf
	zip $(TARGET).zip *.sty *.bib *.bst *.cls Makefile $(TARGET).tex `for i in \`cat $(TARGET).tex |  grep includegra | grep -v ^% | grep -v "^ %" | cut -f 2 -d'{' | cut -f 1 -d'}'\`; do ls $$i ; ls $$i.pdf ; ls $$i-eps-converted-to.pdf ; ls $$i.png ; ls $$i.jpg; done`
	zip $(TARGET).zip `for i in \`cat $(TARGET).tex |  grep animategra | grep -v ^% | grep -v "^ %" | cut -f 3 -d'{' | cut -f 1 -d'}'\`; do for j in \`ls $$i\` ; do ls $$i$$j ; done ; done`
#	zip $(TARGET).zip $(TARGET).pdf *.sty *.bib *.bst *.cls Makefile $(TARGET).tex `for i in \`cat $(TARGET).tex |  grep includegra | grep -v ^% | grep -v "^ %" | cut -f 2 -d'{' | cut -f 1 -d'}'\`; do ls $$i ; done`


$(TARGET).pdf: $(TARGET).tex
	pdflatex $(TARGET).tex

clean:
	-rm $(TARGET).aux $(TARGET).log $(TARGET).nav $(TARGET).out $(TARGET).snm $(TARGET).toc $(TARGET).vrb

purge: clean
	rm $(TARGET).pdf 
	
