echo \\documentclass{minimal}\
     \\usepackage{graphicx}\
     \\begin{document}\
     \\includegraphics{$1}\
     \\end{document} > ZZZtemporarygarbage.tex
latex ZZZtemporarygarbage.tex
dvipng -v -D 100 -T tight -bg Transparent -o ZZZtemporarygarbage.png ZZZtemporarygarbage.dvi
display ZZZtemporarygarbage.png &
#rm ZZZtemporarygarbage.*
