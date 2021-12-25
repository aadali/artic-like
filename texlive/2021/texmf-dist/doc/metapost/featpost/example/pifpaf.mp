% pifpaf.mp
% L. Nobre G.
% 2010

prologues := 1;

beginfig(1);
  numeric u, siz, lena, lenb, lenc, anga, angb, angc, astep;
  u = 7mm;
  siz = 12u;
  lena = 9u;
  lenb = 8u;
  lenc = 4u;
  astep = 1.7;
  z1 = (-lenb,siz-lena);
  z2 = (lenb,siz-lena);
  path circa, circb;
  numeric auxa, auxb, auxc, lend;
  auxa = 2*(lena+2*lenc);
  auxb = 2*lena;
  circa = fullcircle scaled auxa shifted z1;
  circb = fullcircle scaled auxb shifted z2;
  z3 = circa intersectionpoint circb;
  anga = angle(z3-z1);
  angb = 180-angle(z3-z2);
  auxc = lena +-+ (lenb-lenc);
  angc = angle(lenb-lenc,auxc);
  show (anga,angc,angb);
  z4 = z2+(lenc-lenb,auxc);
  z5 = z1+(lenb-lenc,auxc);
  draw z1--z2--z4--z5--cycle;
  draw z1--z2--z3--cycle;
  z6 = z1+lena*dir(anga);
  z7 = z2+lena*dir(180-angb);
  z8 = 0.5[z6,z7];
  lend = siz-lena+auxc;
  z9 = z8+lend*dir(anga-90);
  draw z8--z9;
  draw origin withpen pencircle scaled 3pt;
  numeric i, diamc;
  pair parp, intp, pa, pb;
  diamc = 4*lenc;
  for i=anga step astep until angc:
    parp := z1+lena*dir(i);
    circa := fullcircle scaled diamc shifted parp;
    intp := circa intersectionpoint circb;
    pa := 0.5[parp,intp];
    pb := lend*unitvector((ypart (intp-parp),-xpart (intp-parp)));
    draw pa--(pa+pb);
    draw (-xpart (pa+pb),ypart (pa+pb)) withpen pencircle scaled 2pt;
    draw (pa+pb) withpen pencircle scaled 2pt;
  endfor;
  
endfig;

end.