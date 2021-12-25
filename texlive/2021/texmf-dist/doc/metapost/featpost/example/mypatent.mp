% mypatent.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

verbatimtex 
\documentclass[12pt]{article}
\usepackage{courier}
\begin{document}
{\ttfamily
etex

vardef P@# = rp(V@#) enddef;

picture pecaa, pecab, pecabe;

f := 2*(14,5,10);
Spread := 36;
ahangle := 30;
ahlength := 1.8mm;

pen pendot, penedge, pendist;
penedge= pencircle scaled 1pt;
pendot = pencircle scaled 1.7pt;
pendist= pencircle scaled 0.4pt;

numeric u, marg, leves, elast;
u = 0.1;
marg = 16u;
%leves = 410;
leves = 22u;
elast = 1u;

vardef elastf( expr Ca, Cb ) = Cb + N( Cb-Ca )*elast enddef;

V1 := (0,0,0);
V2 := (0,0,4u);
V3 := (0,16u,4u);
V4 := (4u,16u,4u);
V5 := (4u,16u,8u);
V6 := (0,16u,8u);
V7 := (0,24u,16u);
V8 := (4u,24u,16u);
V9 := (4u,35u,4u);
V10:= (4u,35u,0);
V11:= (16u,23u,0);
V12:= (16u,23u,4u);
V13:= (16u,0,4u);
V14:= (16u,0,0);
V15:= (0,43u,8u);
V16:= (4u,43u,8u);
V17:= (8u,8u,0);
V18:= (4u,43u,16u);
V19:= (0,43u,16u);

beginfig(1);
  path pat[];
  pat1 = P3--P2--P13--P14--P11--P10--P9--P12--P13;
  pat2 = P11--P12;
  pat3 = P9--P4--P5--P6--P7--P8--P5;
  numeric i;
  pat8 = P1-- for i=-90 step 5 until 85:
    rp(V15+(0,8u*cosd(i),8u*sind(i)))..
    endfor P19--P7;
  pat5 = P10-- for i=-90 step 5 until 85:
    rp(V16+(0,8u*cosd(i),8u*sind(i)))..
    endfor P18--P8;
  pat6 = pat1 cutbefore pat3;
  pair auxp,auxP;
  auxp = point 0 of spatialhalfcircle(V15,red,8.05u,true);
  auxP = point 0 of spatialhalfcircle(V16,red,8.05u,true);
  pat9 = rp(V15)--auxp;  
  pat4 = pat8 cutbefore pat9;
  pat7 = auxp--auxP;
  drawoptions( withpen penedge );
  for i=2 upto 7:
    draw pat[i];
  endfor;
  fakehole(V15,red*4u, 2u );
  fakehole(V17,blue*4u, 2u );
  pecab = currentpicture;
endfig;

beginfig(2);
  numeric i;
  for i=1 upto 19:
    V[i+100] := (X(V[i]),-Y(V[i]),Z(V[i]))-(0,marg,0);
  endfor;
  path pat[];
  pat1 = P109--P112--P113--P102--P101--P114--P111--P110--P109--P104--P103
    --P106--P107--P108--P105--P104;
  pat2 = P105--P106;
  pat3 = P113--P114;
  pat4 = P111--P112;
  pat5 = P102--P103;
  pat10 = P101-- for i=-90 step 5 until 85:
    rp(V115+(0,-8u*cosd(i),8u*sind(i)))..
    endfor P119--P107;
  pat6 = P110-- for i=-90 step 5 until 85:
    rp(V116+(0,-8u*cosd(i),8u*sind(i)))..
    endfor P118--P108;
  pair auxp,auxP;
  auxp = point 0 of spatialhalfcircle(V115,red,8.05u,false);
  auxP = point 0 of spatialhalfcircle(V116,red,8.05u,false);
  pat9 = P115--auxp;  
  pat8 = pat10 cutbefore pat9;
  pat7 = auxp--auxP;
  drawoptions( withpen penedge );
  for i=1 upto 8:
    draw pat[i];
  endfor;
  fakehole(V115,red*4u, 2u );
  fakehole(V117,blue*4u, 2u );
  pecabe = currentpicture;
%  draw pecab;
endfig;

beginfig(3);
  V21:= (24u+marg,0,0);
  V22:= (X(V21),Y(V115)+16u,0);
  V23:= (X(V21),Y(V15)-16u,0);
  V24:= (X(V21),Y(V115),0);
  V25:= (X(V21),Y(V15),0);
  V26:= (X(V21),Y(V115),4u);
  V27:= (X(V21),Y(V15),4u);
  V28:= (X(V21),Y(V115)+16u,4u);
  V29:= (X(V21),Y(V15)-16u,4u);
  path pat[];
  pat1 = for i=-90 step 2 until 90:
    rp(V22+(2u*sind(i),-2u*cosd(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V23+(2u*sind(-i),2u*cosd(i),0))...
    endfor
    cycle;
  pat2 = for i=-90 step 2 until 90:
    rp(V28+(2u*sind(i),-2u*cosd(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V29+(2u*sind(-i),2u*cosd(i),0))...
    endfor
    cycle;
  pat3 = for i=-90 step 2 until 90:
    rp(V24+(8u*sind(i),-8u*cosd(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V25+(8u*sind(-i),8u*cosd(i),0))...
    endfor
    cycle;
  pat4 = for i=-90 step 2 until 90:
    rp(V26+(8u*sind(i),-8u*cosd(i),0))...
    endfor
    for i=-90 step 2 until 90:
    rp(V27+(8u*sind(-i),8u*cosd(i),0))...
    endfor
    cycle;
  pair auxp, auxP, auxm, auxM;
  auxp = point 0 of spatialhalfcircle(V24,blue,8.05u,true);
  auxP = point 0 of spatialhalfcircle(V26,blue,8.05u,true);
  auxm = point 0 of spatialhalfcircle(V25,blue,8.05u,false);
  auxM = point 0 of spatialhalfcircle(V27,blue,8.05u,false);
  pat5 = P24--auxp;
  pat6 = P25--auxm;
  pat7 = pat3 cutbefore pat5 cutafter pat6;
  pat8 = auxp--auxP;
  pat9 = auxm--auxM;
  drawoptions( withpen penedge );
%  draw pat1;
  draw pat2;
  draw pat7;
  draw pat4;
  draw pat8;
  draw pat9;
  fakehole(V24,blue*4u, 2u );
  fakehole(V25,blue*4u, 2u );
  pecaa = currentpicture;
%  draw pecab;
%  draw pecabe;
endfig;

beginfig(4);
  V31:= V26-(5u,0,0);
  V32:= 0.25*(V105+V106+V107+V108);
  V33:= 0.25*(V5+V6+V7+V8);
  pair paa, pab, pabe, angs, pn[];
  numeric anxy, angz;
  color veca, vecb, auxc, vecc;
  angs = getanglepair( f-viewcentr );
  anxy = xpart angs;
  angz = ypart angs;
  vecb = eulerrotation( anxy, angz+85.5, 0, red );
  veca = eulerrotation( anxy, angz+83.5, 0, red );
  vecc = eulerrotation( anxy, angz+93, 0, red );
  drawoptions( withpen pendist );
  pn2 = rp( lineintersectplan( V26, blue, f, vecb ) );
  draw P26--pn2;
  label.top( btex 3 etex, pn2 );
  pn3 = rp( lineintersectplan( V27, blue, f, vecc ) );
  draw P27--pn3;
  label.top( btex 3' etex, pn3 );
  auxc = 0.65[V26,V27];
  pn4 = rp( lineintersectplan( auxc, blue, f, vecb ) );
  draw rp(auxc)--pn4;
  label.top( btex 4 etex, pn4 );
  pn5 = rp( lineintersectplan( V116, blue, f, vecb ) );
  draw P116--pn5;
  label.top( btex 5 etex, pn5 );
  pn6 = rp( lineintersectplan( V117+(0,0,4u), blue, f, vecb ) );
  draw rp(V117+(0,0,4u))--pn6;
  label.top( btex 6 etex, pn6 );
  pn7 = rp( lineintersectplan( 0.25*(V103+V104+V105+V106), blue, f, vecb ) );
  draw rp(0.25*(V103+V104+V105+V106))--pn7;
  label.top( btex 7 etex, pn7 );
  pn8 = rp( lineintersectplan( 0.5*(V102+V117+(0,0,4u)), blue, f, vecb ) );
  draw rp(0.5*(V102+V117+(0,0,4u)))--pn8;
  label.top( btex 8 etex, pn8 );
  
  
%   paa = rp((X(V31),Y(V31),leves));
%   pab = rp((X(V33),Y(V33),leves));
%   pabe= rp((X(V32),Y(V32),leves));
  paa = rp( lineintersectplan( V31, blue, f, veca ) );
  pab = rp( lineintersectplan( V33, blue, f, veca ) );
  pabe= rp( lineintersectplan( V32, blue, f, veca ) );
  drawoptions( withpen pendot );
  draw P31;
  draw P32;
  draw P33;
  draw paa;
  draw pab;
  draw pabe;
  drawoptions( withpen pendist );
  draw P31--paa;
  draw P32--pabe;
  draw P33--pab;
  label.top( btex 1 etex, paa );
  label.top( btex 2 etex, pabe );
  label.top( btex 2' etex, pab );
%   draw pecaa;
%   draw pecab;
%   draw pecabe;
  V34:= V26+(marg,0,0);
  V35:= V27+(marg,0,0);
  V36:= elastf(V26,V34 );
  V37:= elastf( V27,V35 );
  V38:= V27+(8u,0,0);
  V39:= V25+(8u,0,0);
  V40:= V27-(8u,0,0);
  V41:= V38+(0,marg,0);
  V42:= V39+(0,marg,0);
  V43:= V40+(0,marg,0);
  V44:= elastf( V38,V41 );
  V45:= elastf( V39,V42 );
  V46:= elastf( V40,V43 );
  V47:= 0.5[V34,V35];
  V48:= 0.5[V41,V42];
  V49:= 0.5[V41,V43];
  draw P26--P36;
  draw P27--P37;
  draw P38--P44;
  draw P39--P45;
  draw P40--P46;
  drawdblarrow P34--P35;
  drawdblarrow P41--P42;
  drawdblarrow P41--P43;
  label.bot( btex C etex, P47 );
  label.rt( btex E etex, P48 );
  label.rt( btex L etex, P49 );
  V50:= V109 - (0,4u,0);
  V51:= V115 + (4u,0,0);
  V52:= V51 + (0,4u,0);
  V53:= elastf( V109,V50 );
  V54:= elastf( V51,V52 );
  V55:= (8u,8u,4u);
  V56:= (8u,8u,8u);
  V57:= elastf( V55,V56 );
  V58:= V15+(4u,0,0);
  V59:= V58+(4u,0,0);
  V60:= elastf( V58,V59 );
  
  
  drawoptions( withpen penedge );
  draw pecaa;
  draw pecab;
  draw pecabe;
  drawoptions( withpen pendist );
  draw P109--P53;
  draw P51--P54;
  drawdblarrow P50--P52;
  label.rt( btex E etex, 0.5[P52,P50] );
  draw P55--P57;
  draw P58--P60;
  drawdblarrow P56--P59;
  label.top( btex D etex, 0.5[P56,P59] );

endfig;
  
verbatimtex 
}
\end{document} 
etex

end;















