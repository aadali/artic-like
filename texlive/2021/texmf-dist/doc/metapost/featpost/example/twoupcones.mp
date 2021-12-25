% twoupcones.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

f := 3*(2,1,0.5);

beginfig(1);
  numeric raya, rayb, heia, heib, rang, marg, u;
  color cena, cenb, vera, verb, auxa, auxb, pora, porb, poha, pohb;
  u = 0.38;
  raya = 1u;
  rayb = 2u;
  heia = 0.65u;
  marg = 0.2u;
  rang = 60;
  heib = heia*rayb/raya;
  cena = (raya*(1+marg),0,0);
  cenb = (-rayb*(1+marg),0,0);
  auxa = cena + raya*(cosd(rang),sind(rang),0);
  auxb = cenb + rayb*(cosd(rang),sind(rang),0);
  pora = 0.5[auxa,cena];
  porb = 0.5[auxb,cenb];
  vera = cena+heia*blue;
  verb = cenb+heib*blue;
  poha = 0.5[cena,vera];
  pohb = 0.5[cenb,verb];
  verygoodcone(true,cenb,blue,rayb,verb);
  draw rp(verb)--rp(cenb)--rp(auxb) dashed evenly;
  label.urt(btex $h_1$ etex, rp(pohb));
  label.urt(btex $R_1$ etex, rp(porb));
  verygoodcone(true,cena,blue,raya,vera);
  draw rp(vera)--rp(cena)--rp(auxa) dashed evenly;
  label.urt(btex $h_2$ etex, rp(poha));
  label.urt(btex $R_2$ etex, rp(pora));
endfig;

verbatimtex 
\end{document} 
etex

end;


