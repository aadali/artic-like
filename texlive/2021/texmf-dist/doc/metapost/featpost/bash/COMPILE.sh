#! /bin/bash
rm -fR mps png
mkdir png
mkdir mps
cd example
mpost --tex=latex graphs.mp
mpost --tex=latex hexagonalgrid.mp
mv graphs.? graphs.??? hexagonalgrid.? hexagonalgrid.??? ../
cd ..
echo \\documentclass[twocolumn]{article}\
     \\usepackage[dvips]{epsfig}\
     \\begin{document}\
     \\title{Exemplifier}\
     \\author{L. Nobre G.}\
     \\date{\\today}\
     \\maketitle\
     \\begin{center}\
     > Exemplifier.tex
for a in example/*.mp
do
  b=`basename $a .mp`
  mpost1.504 --tex=latex example/$b
  for c in $b.?
  do
    echo \\vfill\
      \\begin\{tabular\}\{c\}\
      \\epsfig\{file=$c\,width=0.9\\columnwidth\} \\\\\
      \\texttt\{$c\}\
      \\end\{tabular\} \\\\ >> Exemplifier.tex
  done
done
echo \\end{center}\\end{document} >> Exemplifier.tex
latex Exemplifier.tex
dvips -o Exemplifier.ps Exemplifier.dvi
ps2pdf Exemplifier.ps
#bzip2 Exemplifier.ps
rm Exemplifier.aux
rm Exemplifier.dvi
rm Exemplifier.tex
#mv Exemplifier.ps.bz2 doc/
mv Exemplifier.pdf doc/
cp *.? mps/
echo "<HTML><HEAD><TITLE>MetaPost Examples</TITLE>"\
     "<META NAME=\"KEYWORDS\" content=\"MetaPost, Examples, examples, code\">"\
     "</HEAD><BODY TEXT=\"#000000\" LINK=\"#0000FF\" VLINK=\"#AA00AA\""\
     "ALINK=\"#AA0000\" BGCOLOR=\"#b0c4de\">"\
     "<H1><CODE>MetaPost</CODE> Examples</H1>" > MPexamples.html
for a in example/*.mp
do
  b=`basename $a .mp`
  echo "<P align=\"center\"><a href=\"./example/"$b.mp"\">"\
       $b.mp"</a>" >> MPexamples.html
  for c in $b.?
  do
###############################################################################
    ./laproof.sh $c 
    ./pngmaker.sh $c.eps
###############################################################################
    echo "<img align=\"middle\" src=\"./png/"$c".png\" alt=\""\
         $c "\">" >> MPexamples.html
  done
  echo "</P>" >> MPexamples.html
done
echo "<HR SIZE=4><BR><ADDRESS><A HREF=\"lnobreg at gmail\">L. Nobre G.</A>"\
     "</ADDRESS></BODY></HTML>" >> MPexamples.html
rm *.? 
rm *.?.aux *.?.latex *.?.dvi *.log *.?.mgk *.?.orig *.?.eps
rm example/*.mpx *.mpx
cd doc
pdflatex featpostmanual.tex
pdflatex featpostmanual.tex
#rm -fR ../mps
