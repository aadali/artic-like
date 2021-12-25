echo \\documentclass{minimal}\
     \\usepackage{graphicx}\
     \\DeclareGraphicsRule{*}{mps}{*}{}\
     \\begin{document}\
     \\includegraphics{$1}\
     \\end{document} > ZZZtemporarygarbage.tex
pdflatex ZZZtemporarygarbage.tex
xpdf ZZZtemporarygarbage.pdf &
#rm ZZZtemporarygarbage.*