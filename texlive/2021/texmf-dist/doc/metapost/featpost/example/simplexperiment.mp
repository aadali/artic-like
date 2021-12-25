% simplexperiment.mp 
% L. Nobre G.
% 2007

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{palatino} 
\usepackage{mathpazo}  
\begin{document} 
etex

beginfig(1);
  Spread := 45;
  numeric barh, barl, bart, bardist, bardept, sphr, numlats, i, ie, is; 
  numeric cicon, alpha, fhd, conr, conh, linel;
  color centerb, linpa, linpb, linpc;
  path tmp, lina, cicloi;
  pen thin, thick, aver;
  sphr = 1;
  conr = 0.3;
  barh = 1;
  barl = 7;
  bart = 0.2;
  numlats = 2;
  linel = 5;
  ie = 630;
  is = 5;
  alpha = 40;
  fhd = 9;
  f := 3*(fhd*cosd(alpha)-sphr*sind(alpha),fhd*sind(alpha)+sphr*cosd(alpha),2);
  bardist = sphr*sind(90/numlats);
  bardept = sphr*cosd(90/numlats);
  conh = barh+bardept-sphr;
  cicon = bardept*3.14159/180.0;
  thin = pencircle scaled 0.5pt;
  thick= pencircle scaled 1.5pt;
  aver = pencircle scaled 1.0pt;
  ahlength:=2mm;
  ahangle:=30;
%  kindofcube( WithDash, IsVertex, RefP, AngA, AngB, AngC, LenA, LenB, LenC )
  kindofcube(false,true,(-bardist-bart,0,-bardept),0,0,180,bart,barl,barh );
%  tropicalglobe( expr NumLats, TheCenter, Radius, AngulMom );
  tropicalglobe( numlats, black, sphr, red );
  cicloi = rp((0,0,sphr)) for i=is step is until ie:
    ..rp((0,-sphr*sind(i)-cicon*i,sphr*cosd(i))) endfor;
  draw cicloi;
  kindofcube( false, true, (bardist,0,-bardept), 0,0,180, bart, barl, barh );
  centerb = 2*sphr*(-sind(alpha),cosd(alpha),0);
  linpb = -blue*(bardept+barh);
  linpa = linpb+green*linel;
  linpc = linpb+linel*(green*cosd(alpha)-red*sind(alpha));
  draw rp(linpa)--rp(linpb)--rp(linpc);
  angline( linpa, linpc, linpb, 0.7*linel, btex $\alpha$ etex scaled 3)(rt);
  verygoodcone(false,centerb-(sphr+conh)*blue,blue,conr,centerb-sphr*blue);
  tmp = rigorousfearpath( centerb, sphr );
  unfill tmp;
  draw tmp;
  drawoptions( dashed evenly );
  draw goodcirclepath( black, blue, sphr );
  draw goodcirclepath( centerb, blue, sphr );
  


endfig;

verbatimtex 
\end{document} 
etex

end;
