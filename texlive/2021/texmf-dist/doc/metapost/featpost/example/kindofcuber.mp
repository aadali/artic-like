% kindofcube.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

verbatimtex \documentclass{article}\begin{document}{\sffamily etex

beginfig(1);
  f := (4,1.2,0.5);
  Spread := 120;
  numeric lena, lenb, lenc, aa, ab, ac, angray, micr, labscal;
  pair zeror;
  color axa, axb, axc, returnedcorner, dummyc, auxx, auxy;
  lena=0.3;
  lenb=0.6;
  lenc=0.9;
  aa=130;
  ab=32;
  ac=67;
  labscal=1.6;
  angray=0.85;
  micr=0.01;
  zeror=rp(black);
  cartaxes(1.2,1,1.15);
  axa = ( cosd(aa), sind(aa), 0 );
  axb = axa*cosd(ab)+blue*sind(ab);
  auxx= ( cosd(aa+90), sind(aa+90), 0);
  auxy = ncrossprod( axb, auxx );
  axc = auxx*cosd(ac)+auxy*sind(ac);
  draw zeror--rp(axa);
  draw zeror--rp(axb);
  draw zeror--rp(-0.6*axa);
  draw zeror--rp(axc);
  draw zeror--rp(auxx);
  label.rt( btex $a$ etex scaled labscal, rp(axa) );
  label.rt( btex $b$ etex scaled labscal, rp(axb) );
  label.lft( btex $c$ etex scaled labscal, rp(auxx) );
  label.lft( btex $d$ etex scaled labscal, rp(axc) );
  angline( axa,red,black,angray,btex $\alpha_1$ etex scaled labscal,bot );
  angline( axb,axa,black,angray,btex $\alpha_2$ etex scaled labscal,rt );
  angline( axc,auxx,black,angray,btex $\alpha_3$ etex scaled labscal,lft );
  squareangline( -axa,auxx,black,0.2);
  drawoptions( withpen pencircle scaled 1.7pt withcolor 0.65blue );
  kindofcube( false, true, black, aa, ab, ac, lena, lenb, lenc );
  drawoptions();
  dotlabel.bot( btex $\vec{o}$ etex scaled labscal, zeror );
  label.bot( btex $l_1$ etex scaled labscal, rp(lena*axb) );
  label.bot( btex $l_2$ etex scaled labscal, rp(0.5*lenb*axc) );
  label.ulft( btex $l_3$ etex scaled labscal,rp(0.5*lenc*ncrossprod(axb,axc)));
endfig;
  
verbatimtex 
\end{document} 
etex

end.

  
