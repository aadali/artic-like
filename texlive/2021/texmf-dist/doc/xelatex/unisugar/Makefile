PACKAGE=unisugar
SOURCES=${PACKAGE}.tex ${PACKAGE}.sty *.png Makefile README

${PACKAGE}.pdf: ${SOURCES}
	yes q | xelatex ${PACKAGE}.tex >/dev/null
	yes q | xelatex ${PACKAGE}.tex >/dev/null
	yes q | xelatex ${PACKAGE}.tex >/dev/null

wrap: ${PACKAGE}.zip 

${PACKAGE}.zip: ${PACKAGE}.pdf $(SOURCES) *.png
	zip $@ $?

clean:
	rm -f *.aux *.toc *~ *.log *.out *.pdf
