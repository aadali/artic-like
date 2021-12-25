#!/bin/bash 
echo \\documentclass{minimal} > $1.latex 
echo \\usepackage[dvips]{epsfig} >> $1.latex 
echo \\begin{document} >> $1.latex 
echo \\epsfig{file=$1} >> $1.latex 
echo \\end{document} >> $1.latex 
latex $1.latex
dvips -E -o $1.eps $1.dvi