% rungekuttasecond.mp
% L. Nobre G.
% 2008

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

Spread := 30;
f := (10,0,0);

def vecfunc( expr posit )=
  ( (0,1,Z(posit)-Y(posit)) )
enddef;
  
beginfig(2);
%  cartaxes( 2.4, 2.4, 2.4 );
  numeric grids, size, i, j, ijmax, diffstep, len, h;
  color actpos, one, two, aux, graycol, stap, midp, finp ;
  pair yaxp, oaxp, xaxp;
  graycol = 0.67white;
  h = 3;
  grids = 0.15;
  size = 0.14;
  diffstep = 0.0125;
%  len = 320;
  len = 240;
  ijmax = floor( 2.0/grids );
  for j=-ijmax upto ijmax:
    for i=-ijmax upto ijmax:
      actpos := grids*(0,i,j);
      aux := 0.5*size*N( vecfunc( actpos ) );
      one := actpos-aux;
      two := actpos+aux;
      draw rp(one)--rp(two) withcolor graycol;
    endfor;
  endfor;
%   stap = (0,-2,-1.05);
%   draw fieldlinepath( len, stap, diffstep, vecfunc ) withcolor blue;
  stap = (0,-1.5,-0.65);
  draw fieldlinepath( len, stap, diffstep, vecfunc ) withcolor blue;
  midp = stap + 0.5*h*vecfunc( stap );
  finp = stap + h*vecfunc( midp );
  draw rp(stap)--rp(midp);
  draw rp(stap)--rp(finp);
%   drawoptions( withcolor blue );
%   draw fieldlinepath( len, (0,-2,-1.25), diffstep, vecfunc );
%   draw fieldlinepath( len, (0,-2,-1.2), diffstep, vecfunc );
%   draw fieldlinepath( len, (0,-2,-1.15), diffstep, vecfunc );
%   draw fieldlinepath( len, (0,-2,-1.1), diffstep, vecfunc );
%   draw fieldlinepath( len, (0,-2,-1.05), diffstep, vecfunc );
  yaxp = rp((0,-2-grids,2+grids));
  oaxp = rp((0,-2-grids,-2-grids));
  xaxp = rp((0,2+grids,-2-grids));
  drawdblarrow yaxp--oaxp--xaxp;
  label.bot(btex $x$ etex, xaxp );
  label.lft(btex $y$ etex, yaxp );
  dotlabel.bot(btex $x_0$ etex, rp((0,Y(stap),-2-grids)) );
  dotlabel.lft(btex $y_0$ etex, rp((0,-2-grids,Z(stap))) );
  dotlabel.bot(btex $x_0+\frac{h}{2}$ etex, rp((0,Y(midp),-2-grids)) );
  dotlabel.lft(btex $y_0+\frac{h}{2}f_0$ etex, rp((0,-2-grids,Z(midp))) );
  dotlabel.bot(btex $x_0+h$ etex, rp((0,Y(finp),-2-grids)) );
  dotlabel.lft(btex $y_0+hf_{1/2}$ etex, rp((0,-2-grids,Z(finp))) );
endfig;

numeric u, grids, i, j, thex, they, xmax, ymax, vertn, horin, rescal, size;
pair actpos, direct, one, two;

u := 0.5mm;
vertn := 20;
horin := 41;
size := 3u; 
grids := 8u;
ymax := 0.5*(vertn-1);
xmax := 0.5*(horin-1);
rescal := 0.15;

beginfig(1);
  pickup pencircle scaled u;
  for i=-xmax upto xmax:	
    for j=-ymax upto ymax:
      actpos:=grids*(i,j);
      thex:=i*rescal;
      they:=j*rescal;
      direct:=unitvector((1,they-thex));
      one:=actpos-size*direct;
      two:=actpos+size*direct;
      draw one--two;
    endfor; 
  endfor;
endfig;

verbatimtex 
\end{document} 
etex

end;
