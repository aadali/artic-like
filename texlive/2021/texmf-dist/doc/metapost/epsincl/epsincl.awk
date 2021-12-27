# made in BOP s.c., bop@bop.com.pl; public domain software
# history: 20.08.1999: ver. 0.1, pre-release
#          05.01.2000: ver. 0.2, first release
BEGIN {
  hd["0"]=0; hd["1"]=1; hd["2"]=2; hd["3"]=3; hd["4"]=4
  hd["5"]=5; hd["6"]=6; hd["7"]=7; hd["8"]=8; hd["9"]=9
  hd["a"]=hd["A"]=10; hd["b"]=hd["B"]=11; hd["c"]=hd["C"]=12
  hd["d"]=hd["D"]=13; hd["e"]=hd["E"]=14; hd["f"]=hd["F"]=15
  #
  hb["0"]="0000"; hb["1"]="0001"; hb["2"]="0010"; hb["3"]="0011"
  hb["4"]="0100"; hb["5"]="0101"; hb["6"]="0110"; hb["7"]="0111"
  hb["8"]="1000"; hb["9"]="1001"; hb["a"]="1010"; hb["b"]="1011"
  hb["c"]="1100"; hb["d"]="1101"; hb["e"]="1110"; hb["f"]="1111"
  #
  for (i in hb) bh[hb[i]]=i
  #
  hb["A"]="1010"; hb["B"]="1011"; hb["C"]="1100"
  hb["D"]="1101"; hb["E"]="1110"; hb["F"]="1111"
  #
  comments=1
}

/^[^%]/ && (comments==1) {comments=0; flush_prolog()}

comments==1 {process_prolog()}

(l==7) {l=0}

(l==6) {
  if ($0~/ *closepath/) {l=7; IncludeEPS(Anum)} else l=-6
}

(l==5) {
  if (is_lineto()) {l=6}
  else if (is_curveto()) {l=6}
  else if ($0~/ *closepath/) {l=7; IncludeEPS(Anum)}
  else l=-5
}

(l==4) {
  if (is_lineto()) {l=5; Aulx=$1; Auly=$2}
  else if (is_curveto()) {l=5; Aulx=$5; Auly=$6}
  else l=-4
}

(l==3) {
  if (is_lineto()) {l=4; Aurx=$1; Aury=$2}
  else if (is_curveto()) {l=4; Aurx=$5; Aury=$6}
  else l=-3
}

(l==2) {
  if (is_lineto()) {l=3; Adrx=$1; Adry=$2}
  else if (is_curveto()) {l=3; Adrx=$5; Adry=$6}
  else l=-2
}

(l==1) {
  if ($0~/^newpath -?[0-9\.]+ -?[0-9\.]+ moveto$/) {
    l=2
    Adlx=$2; Adly=$3
  } else l=-1
}

(l<0) {
  print "Something strange in EPS colour (" (-l) ")" > "CON"
  print NR, is_lineto(), is_curveto() > "CON"
  print $0 > "CON"
  exit
  l=0
}

/setrgbcolor$/ && (n>0) && (comments==0) {
  for (i in Ename)
    if ((Ecolr[i]==$1) && (Ecolg[i]==$2) && (Ecolb[i]==$3)) {
      Anum=i; l=1
    }
}

(l==0) && (comments==0) {gsub(/fshow/, "nshow"); print}

function hexval(x) {return hd[substr(x,1,1)]*16+hd[substr(x,2,1)]}

function process_prolog(){
  if (/^%\*Font:/) {
    f_ch[$2 "\001" $4]=fix_chars(f_ch[$2 "\001" $4],$5)
  } else
  if (/^%EPS/) {
    Ecolr[++n]=$2+0; Ecolg[n]=$3+0; Ecolb[n]=$4+0
    Ename[n]=$5 ""
    Ebbl[n]=$6+0; Ebbd[n]=$7+0; Ebbr[n]=$8+0; Ebbu[n]=$9+0
  } else
    prolog[++prolog[0]]=$0
}

function flush_prolog( j,l){
  for (l=1; l<=prolog[0]; ++l) {
   if (prolog[l]~/^%%EndProlog/) j=flush_fonts();
   print prolog[l]
  }
  if (j>0) print "/nshow {gsave 10 div dup scale 10 fshow grestore} def"
  if (n>0) print "/MEPSDict 50 dict def"
}

function flush_fonts( i,j) {
  j=0
  for (i in f_ch) {
    j++
    split(i,a,"\001")
    printf("%%*Font: %s 10.0 %s %s\n", a[1], a[2], out_chars(f_ch[i]))
  }
  return j
}

function fix_chars(f,c,  a) {
 split(c,a,":")
 return mix_chars((f=="" ? sprintf("%0256d", 0) : f), a[2], hexval(a[1]))
}

function mix_chars(f,c,x,  h,k,l) {
 for (k=1; k<=length(c); k++) {
  h=hb[substr(c,k,1)]
  for (l=1; l<=length(h); l++)
    if (substr(h,l,1)=="1") f=substr(f,1,x+4*k+l-5) "1" substr(f,x+4*k+l-3)
 }
 return f
}

function out_chars(f,  i,j,k,g) {
  while ((i<256) && (substr(f,i,1)=="0")) ++i
  g=""
  for (j=i; j<=length(f); j+=4) {
    k=sprintf("%-4s", substr(f,j,4)); gsub(/ /,"0",k); g=g bh[k]
  }
  sub (/0+$/,"",g)
  return sprintf("%02X", (i-1)) ":" g
}

function is_lineto() {
 return ($0~/^-?[0-9\.]+ -?[0-9\.]+ lineto$/)
}

function is_curveto() {
 return ($0~/^-?[0-9\.]+ -?[0-9\.]+ -?[0-9\.]+ -?[0-9\.]+ -?[0-9\.]+ -?[0-9\.]+ curveto$/)
}

function IncludeEPS(i,  comments, s) {
  print "%BeginInclude"
  print "MEPSDict begin /MEPSSave save def gsave"
  print "count /MEPSocount exch def /MEPSdcount countdictstack def"
  print "/showpage {} def /erasepage {} def /copypage {} def"
  print "/letter {} def /note {} def /a4 {} def /legal {} def"
  printf ("[ %f %f %f %f %f %f ] matrix invertmatrix\n",
    Ebbr[i]-Ebbl[i], 0, 0, Ebbu[i]-Ebbd[i], Ebbl[i], Ebbd[i])
  printf ("[ %f %f %f %f %f %f ] matrix concatmatrix\n",
    Adrx-Adlx, Adry-Adly, Aulx-Adlx, Auly-Adly, Adlx, Adly)
  print "concat"
  print "0 setlinecap 0 setlinejoin 10 setmiterlimit [] 0 setdash"
  print "1 setlinewidth 0 setgray newpath"
  print "%StartInclude"
#
  comments=1
  while (getline s < Ename[i] >0) 
    if ((s!~/^%[%!\*]/) || (comments==0)) {
      comments=0; gsub(/fshow/, "nshow", s); print s
    }
  close(Ename[i])
#
  print "%StopInclude"
  print "count MEPSocount sub {pop} repeat"
  print "countdictstack MEPSdcount sub {end} repeat"
  print "grestore MEPSSave restore end"
  print "%EndInclude"
}


