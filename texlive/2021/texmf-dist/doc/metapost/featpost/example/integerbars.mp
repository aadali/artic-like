% integerbars.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

prologues := 1;

def randcolor =
  (0.1+uniformdeviate(0.8),0.1+uniformdeviate(0.8),0.1+uniformdeviate(0.8))
enddef;

beginfig(1);
  numeric u, ltl, ltr, lsl, lsr, lbmin, lbmax, lbtop, i, j, k, a;
  numeric hist[], hclass, aaux, r;
  color pickc;
  pair tcross;
  pickup pencircle scaled 1.5pt;
  u = 5mm;
  r = 5mm;
  ltl = 9;
  ltr = 9;
  lsl = 4;
  lsr = 1;
  lbmin = 4;
  lbtop = 7;
  lbmax = ltl+ltr-1;
  hclass = 2;
  if lbtop > lbmax: lbtop := lbmax; fi;
  draw origin;
  for a=0 step hclass until 180-hclass:
    hist[a] = 0;
  endfor;
  for i=lbmin upto lbtop:
    pickc := randcolor;
    drawoptions( withcolor pickc );
    draw (i*u,0);
    for j=lsl upto ltl:
      for k=lsr upto ltr:
	if (j+k > i) and (j < i+k) and (k < i+j):
	  tcross := radialcross( origin, j*u, (i*u,0), k*u, true );
	  draw tcross;
	  aaux := angle( tcross );
	  a := 0;
	  forever:
	    a := a + hclass; 
	    exitif aaux < a;
	  endfor;
	  hist[a] := incr( hist[a] );
	fi;
      endfor;
    endfor;
  endfor;
  fill origin-- for a=hclass step hclass until 180-hclass:
    (r*dir(a)*hist[a-hclass])--(r*dir(a)*hist[a])-- endfor cycle;
endfig;

end.