#! /bin/bash
#cd example
#mpost --tex=latex graphs.mp
#mv graphs.? graphs.??? ../
#cd ..
echo "<HTML><HEAD><TITLE>MetaPost Examples</TITLE>"\
     "<META NAME=\"KEYWORDS\" content=\"MetaPost, Examples, examples, code\">"\
     "</HEAD><BODY TEXT=\"#000000\" LINK=\"#0000FF\" VLINK=\"#AA00AA\""\
     "ALINK=\"#AA0000\" BGCOLOR=\"#b0c4de\">"\
     "<H1><CODE>MetaPost</CODE> Examples</H1>" > MPexamples.html
for a in example/*.mp
do
  b=`basename $a .mp`
  #mpost --tex=latex example/$b
  echo "<P align=\"center\"><a href=\"./example/"$b.mp"\">"\
       $b.mp"</a>" >> MPexamples.html
  for c in $b.?
  do
    ./laproof.sh $c 
    ./pngmaker.sh $c.eps
    echo "<img align=\"middle\" src=\"./png/"$c".png\" alt=\""\
         $c "\">" >> MPexamples.html
  done
  echo "</P>" >> MPexamples.html
done
echo "<HR SIZE=4><BR><ADDRESS><A HREF=\"lnobreg at gmail\">L. Nobre G.</A>"\
     "</ADDRESS></BODY></HTML>" >> MPexamples.html
#rm *.? 

