% kindofcube.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

Spread := 30;
f := 5.4*(1.5,0.5,1);

beginfig(1);
  numeric gridstep, sidenumber, i, j, coord, aa, ab, ac;
  color pa;
  gridstep = 0.9;
  sidenumber = 10;
  coord = 0.5*sidenumber*gridstep;
  for i=0 upto sidenumber:
    for j=0 upto sidenumber:
      pa := (-coord+j*gridstep,-coord+i*gridstep,0);
      aa := uniformdeviate(360);
      ab := uniformdeviate(180);
      ac := uniformdeviate(90);
      kindofcube( false, false, pa, aa, ab, ac, 0.4, 0.4, 0.9 );
    endfor;
  endfor;
endfig;

beginfig(2);
  f := (4,1.0,0.5);
  Spread := 120;
  numeric lena, lenb, lenc, aa, ab, ac, angray, micr;
  pair zeror;
  color axa, axb, axc, returnedcorner, dummyc, auxx, auxy;
  lena=0.3;
  lenb=0.6;
  lenc=0.9;
  aa=130;
  ab=20;
  ac=70;
  angray=0.7;
  micr=0.01;
  zeror=rp(black);
  cartaxes(1,1,1);
  axa = ( cosd(aa), sind(aa), 0 );
  axb = axa*cosd(ab)+blue*sind(ab);
  auxx= ( cosd(aa+90), sind(aa+90), 0);
  auxy = ncrossprod( axb, auxx );
  axc = auxx*cosd(ac)+auxy*sind(ac);
  draw zeror--rp(axa);
  draw zeror--rp(axb);
  draw zeror--rp(-axa);
  draw zeror--rp(axc);
  draw zeror--rp(auxx);
  label.rt( btex $a$ etex, rp(axa) );
  label.rt( btex $b$ etex, rp(axb) );
  label.lft( btex $c$ etex, rp(auxx) );
  label.lft( btex $d$ etex, rp(axc) );
  label.bot( btex $o$ etex, zeror );
  angline( axa,red,black,angray,btex $\alpha_1$ etex,bot );
  angline( axb,axa,black,angray,btex $\alpha_2$ etex,rt );
  angline( axc,auxx,black,angray,btex $\alpha_3$ etex,lft );
  squareangline( -axa,auxx,black,0.2);
  kindofcube( true, true, black, aa, ab, ac, lena, lenb, lenc );
  draw zeror--rp(blue) dashed evenly;
  dotlabel.bot( btex $l_1$ etex, rp(lena*axb) );
  dotlabel.rt( btex $l_2$ etex, rp(lenb*axc) );
  dotlabel.top( btex $l_3$ etex, rp(lenc*ncrossprod(axb,axc) ) );
endfig;
  
verbatimtex 
\end{document} 
etex

end.

  
