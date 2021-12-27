MAIN=			drv-guide

$(MAIN).pdf:		$(MAIN).ps
			ps2pdf $(MAIN).ps
			pdfopt $(MAIN).pdf $(MAIN)-opt.pdf
			mv $(MAIN)-opt.pdf $(MAIN).pdf

$(MAIN).ps:		$(MAIN).tex $(MAIN)-delayed.mp
			latex $(MAIN)
			latex $(MAIN)
			latex $(MAIN)
			dvips $(MAIN) -o

$(MAIN)-proof.pdf:	$(MAIN)-proof.ps
			ps2pdf $(MAIN)-proof.ps

$(MAIN)-proof.ps:	$(MAIN)-delayed.mp
			latex $(MAIN)-proof
			dvips $(MAIN)-proof -o

$(MAIN)-delayed.mp:	$(MAIN).mp
			$(RM) $(MAIN)-delayed.mp
			mpost $(MAIN)
			mpost $(MAIN)
			mpost $(MAIN)

.SILENT clean:
			$(RM) drv.mpx
			$(RM) $(MAIN).[0-9]*
			$(RM) $(MAIN).aux
			$(RM) $(MAIN).dvi
			$(RM) $(MAIN).log
			$(RM) $(MAIN).nav
			$(RM) $(MAIN).out
			$(RM) $(MAIN).toc
			$(RM) $(MAIN)-delayed.mp
			$(RM) $(MAIN)-delayed.mpx
			$(RM) $(MAIN)-proof.aux
			$(RM) $(MAIN)-proof.dvi
			$(RM) $(MAIN)-proof.log
			$(RM) $(MAIN)-proof.tex
			$(RM) $(MAIN)-proof.toc

cleanall:		clean
			$(RM) $(MAIN).pdf
			$(RM) $(MAIN).ps
			$(RM) $(MAIN)-proof.pdf
			$(RM) $(MAIN)-proof.ps
