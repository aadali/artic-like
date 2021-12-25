% 38.mp
% J. Dinis
% 2004

def rp(expr p, r)=
	(p + (r * (.5+uniformdeviate 0.5),0) rotated (uniformdeviate 360))
enddef;

vardef eight(expr r) = 
	save	a;
	pair	a;
	a = rp(origin, r);
	(a{dir 15} 
	.. rp((0.6,0.5),r){up} 
	.. rp((0,1),r){left}
	.. rp((-0.6,0.5),r){down} 
	.. rp(origin,r){dir -15} 
	.. rp((0.6,-0.6),r){down} 
	.. rp((0,-1.2),r){left}
	.. rp((-0.6,-0.6),r){up} 
	.. a{dir 15})
enddef;

vardef three(expr r) = 
	save	a;
	path	a;
	a = eight(r);
	((reverse (subpath (0,2.35) of a)) -- (subpath (4,6.5) of a))
enddef;

def scaledto(expr p, wd)=
	p shifted (-center p) scaled (wd / (xpart (urcorner p - llcorner p)))
	shifted (105mm,148.5mm)
enddef;


beginfig(1);
  color	c;
  numeric	n,e;
  pickup pencircle scaled 1pt;
  n := 33;
  for i=n step -2 until 0:
    e := 1 / (i+1);
    pickup pencircle scaled (e * pt);
    for j=1 upto (i+1):
      c := (i/n)[.75blue,white];
      draw eight(i/n) scaled 100 shifted (140,0) withcolor c;
      draw three(i/n) scaled 100 shifted (0,0) withcolor c;
    endfor;
  endfor;
  currentpicture := scaledto(currentpicture, 24cm);
endfig;

end
