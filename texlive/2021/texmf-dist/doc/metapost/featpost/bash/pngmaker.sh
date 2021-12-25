outf=`basename $1 eps`png
prefwidth=160
prefheigh=144
prefratio=$(( 1000 * $prefwidth / $prefheigh ))
tmpf=`mktemp pngmaker.XXXXXX`
epstopdf --nogs -o $tmpf $1
epswidth=`identify $tmpf | cut -d " " -f 3 | cut -d x -f 1`
epsheigh=`identify $tmpf | cut -d " " -f 3 | cut -d x -f 2`
aspectra=$(( 1000 * $epswidth / $epsheigh ))
if test $aspectra -ge $prefratio
then resol=$(( 116 * $prefwidth / $epswidth )) 
else resol=$(( 116 * $prefheigh / $epsheigh ))
fi
echo "                                      Resolution:" $resol "dpi."
gs -q -sDEVICE=pngalpha -r$resol -dNOPAUSE -sOutputFile=$outf -dTextAlphaBits=4 -dGraphicsAlphaBits=4 $tmpf </dev/null
mv $outf png/
rm $tmpf

