% simplependulum.mp
% L. Nobre G.
% 2012

verbatimtex \documentclass{article}\begin{document} etex

prologues := 3;

beginfig(1);
  numeric u, len, hei, ang, maxa, hala, num, signa, marga, margb, ste;
  u = 2.2mm;
  len = 15u;
  hei = 4u;
  maxa = angle( len-hei, len +-+ (len-hei) );
  num = 8;
  hala = 0.7*maxa;
  signa = -1;
  marga = 2u;
  margb = 1u;
  ste = maxa/num;
  z1 = len*(signa*sind(maxa),-cosd(maxa));
  z7 = len*(-signa*sind(hala),-cosd(hala));
  z2 = (0,-len);
  z6 = z1 + (signa*marga,0);
  z5 = z6 + (signa*margb,0);
  z3 = (x6,y2);
  z4 = (x5,y2);
  draw origin--z1;
  draw origin--z2 dashed evenly;
  draw z2--z4 dashed evenly;
  draw z1--z5 dashed evenly;
  path mov, penda, pendb;
  mov = z1 for ang=signa*(maxa-ste) step -signa*ste until -signa*hala:
      ..len*(sind(ang),-cosd(ang)) endfor;
  penda = fullcircle scaled 2u shifted z1;
  pendb = fullcircle scaled 2u shifted z2;
  fill penda withcolor 0.5white;
  draw penda;
  unfill pendb;
  draw pendb dashed evenly;
  drawarrow mov dashed evenly;
  drawdblarrow z3--z6;
  label( btex $h$ etex, 0.5[z4,z5] );
  label( btex $A$ etex, z1+down*(u+margb) );
  label( btex $B$ etex, z2+down*(u+margb) );
  label.top( btex $O$ etex, origin );  
endfig;

verbatimtex \end{document} etex

end.
