% naotocanacalote.mp
% L. Nobre G.
% 2012

verbatimtex \documentclass{article}\begin{document} etex
  
beginfig(1);
  numeric N, h, R, m, i, l, d, X, angl;
  N = 2;
  R = 13mm;
  h = N*R;
  m = 8;
  l = 2R;
  d = 4mm;
  angl = 50;
  path calotes, queda;
  if N=1:
    calotes = (0,h)... for i=1 upto m-1:
	(R*(sind(i*90/m),cosd(i*90/m)))...
    endfor {down}(R,0);
    X = sqrt(2)*R;
    queda = (-R,h+d)---(0,h+d){right}...{unitvector(1,-sqrt(2))}(X,d);
  elseif N=2:
    calotes = (0,h)... for i=1 upto 2*m-1:
	(R*(sind(i*90/m),1+cosd(i*90/m)))...
    endfor {left}origin;
    X = 2R;
    queda = (-R,h+d)---(0,h+d){right}...{unitvector(1,-2)}(X,d);
  else:
    calotes = (0,h)... for i=1 upto 2*m-1:
	(R*(sind(i*90/m),2+cosd(i*90/m)))...
    endfor for i=1 upto m-1:
	(R*(sind(-i*90/m),cosd(i*90/m)))...
    endfor {down}(-R,0);
    X = sqrt(6)*R;
    queda = (-R,h+d)---(0,h+d){right}...{unitvector(1,-sqrt(6))}(X,d);
  fi;
  drawoptions( withpen pencircle scaled (0.75mm) );
  draw (-l,h)--calotes;
  draw (-R,0)--(X+R,0);
  drawoptions();
  draw fullcircle scaled (2d) shifted (0,h+d);
  draw fullcircle scaled (2d) shifted (X,d);
  drawarrow (0,h+d)--(R,h+d);
  label.top( btex $\vec{v}_{\mathrm{min}}$ etex, (0.6R,h+d) );
  drawdblarrow origin--(0,h);
  label.lft( btex $h$ etex, (0,h-0.5R) );
  drawdblarrow (0,d)--(X,d);
  label.top( btex $x$ etex, (0.5X,d) );
  drawarrow (0,h-R)--(R*cosd(angl),h+R*(sind(angl)-1));
  label.ulft( btex $R$ etex, (0.5*R*cosd(angl),h-R+0.5*R*sind(angl)) );
  drawoptions( dashed evenly );
  draw queda;
endfig;

verbatimtex \end{document} etex
  
end.