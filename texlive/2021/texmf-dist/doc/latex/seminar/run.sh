#!/bin/bash
#latex semsamp1 && dvips -T 297mm,210mm semsamp1 && ps2pdf -dALLOWPSTRANSPARENCY semsamp1.ps
latex semsamp1 && dvips semsamp1 && ps2pdf -dALLOWPSTRANSPARENCY semsamp1.ps
latex semsamp2 && dvips -T 297mm,210mm semsamp2 && ps2pdf -dALLOWPSTRANSPARENCY semsamp2.ps
latex semsamp3 && dvips -T 297mm,210mm semsamp3 && ps2pdf -dALLOWPSTRANSPARENCY semsamp3.ps
rm *.log *.aux *.dvi *.ps *.tmp *.out *.los


