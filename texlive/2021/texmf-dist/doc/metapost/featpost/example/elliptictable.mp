% elliptictable.mp
% L. Nobre G. 
% 2006

verbatimtex \documentclass{article} \begin{document} etex

beginfig(1);
  numeric u, l, mar, h, i, s, ai, bi, ao, bo, d;
  u = 0.85mm;
  l = 125u;
  d = 0.5l;
  mar = 12.0u;
  h = 70u;
  s = 5; 
  z1 = (l-2mar,0);
  z2 = (l-mar,0);
  z3 = (l,0);
  z5 = (l-2mar,d-mar);
  z6 = (l-2mar,d);
  z7 = (l-mar,d);
  z8 = (l,d);
  z10= (l-2mar-h,d-mar);
  z11= (l-2mar-h,d);
  z12= (0,d);
  z13= (0,d-mar);
  ai = l-2mar;
  bi = (d-mar)/(1+-+((l-2mar-h)/ai));
  z14= (ai+-+bi,0);
  ao = l-mar;
  bo = d/(1+-+((l-2mar-h)/ao));
  z15= (ao+-+bo,0);
  path ellin, elout, cutpe, cutpb, bordo;
  ellin = z1 for i=s step s until 90: ..(ai*cosd(i),bi*sind(i)) endfor;
  elout = z2 for i=s step s until 90: ..(ao*cosd(i),bo*sind(i)) endfor;
  cutpe = (1.1[z11,z10])--(1.1[z10,z11]);
  cutpb = z1--z6;
  bordo = elout cutafter cutpe cutbefore cutpb;
  z4 = point 0 of bordo;
  z9 = bordo intersectionpoint (z5--z13);
  draw z3--z8--z6--z4;
  draw z5--z9;
  draw z2--z7;
  draw ellin--z12--z11--(reverse bordo)--z1;
  draw z6--z11--z10;
  draw currentpicture xscaled -1;
  pen pendot;
  pendot = pencircle scaled 2u;
  z22= (l,d-mar*ai/bi);
  z23= (l-2mar,d-mar*ai/bi);
%  draw z22{up}..z7{left}..{down}z23;
  z24= z7+down*mar;
  draw z22--z24--z23;
  drawoptions( dashed evenly );
  z17 = (ai,bi)*sqrt(2.0)*0.5;
  z18 = z17 - h*unitvector((bi,ai));
  draw z17--z18;
  z19 = (0,-y18-x18*bi/ai);
  z20 = z19+whatever*(ai,-bi);
  z20 = z17-whatever*(bi,ai);
  z21 = z20-l*unitvector((ai,-bi));
  draw z20--z21;  
  drawoptions();
  draw z14 withpen pendot;
  draw z15 withpen pendot;
  draw z17 withpen pendot;
%  dotlabels.bot( 17, 18, 20, 21 );
  picture tmppic;
  tmppic = currentpicture;
  currentpicture := nullpicture;
  z16= (105mm,145mm);
  draw tmppic rotated 90 shifted z16;
  draw tmppic rotated -90 shifted z16;
  label.rt( decimal((x14)/u), z14 rotated 90 shifted z16 );
  label.lft( decimal((x15)/u), z15 rotated 90 shifted z16 );
  label.rt( decimal((x17)/u), z17 rotated 90 shifted z16 );
  tmppic := currentpicture;
  currentpicture := nullpicture;
  draw tmppic shifted (-z16) rotated (-90) shifted z16;
endfig;

verbatimtex \end{document} etex

end;
