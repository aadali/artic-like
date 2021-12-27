% simplecenterofmass.mp
% L. Nobre G. 
% 2012 

verbatimtex \documentclass{article}\begin{document} etex;

beginfig(0); 						
	numeric u, i, j, s, as, bs, cs, ds;
	pen pena, penb, penc, pend;
	u = 1.5mm;
	s = 7.1u;
	as = 2u;
	bs = 3u;
	cs = 5u;
	ds = 0.8u;
	pena = pencircle scaled as;
	penb = pencircle scaled bs;
	penc = pencircle scaled cs;
	pend = pencircle scaled ds;
	linecap := squared;
	for i=1 upto 4:
	    for j=1 upto 3:
		draw ((i-1)*s,(j-1)*s)--
		     ((i-1)*s,j*s)--
		     (i*s,j*s)--
		     (i*s,(j-1)*s)--cycle withcolor 0.7white;
	    endfor;
	endfor;
	z1 = (s,3s);
	z2 = (3s,3s);
	z3 = (4s,3s);
	z4 = (0,0);
	z5 = (1s,2s);
	z6 = (2s,2s);
	z7 = (3s,2s);
	z8 = (1s,1s);
	z9 = (2s,1s);
	linecap := rounded;
	pickup pena;
	draw z1;
	draw z2;
	pickup penb;
	draw z3;
	pickup penc;
	draw z4;
	pickup pend;
	for i=5 upto 9:
	    draw z[i];
	endfor;
	label.urt(btex 1 kg etex, z1);
	label.urt(btex 1 kg etex, z2);
	label.rt(btex 2 kg etex, z3+(bs/2,0));
	label.urt(btex 8 kg etex, z4+(0,cs/2));
	label.urt(btex A etex, z5);
	label.urt(btex B etex, z6);
	label.urt(btex C etex, z7);
	label.urt(btex D etex, z8);
	label.urt(btex E etex, z9);
endfig;

beginfig(1);
  numeric u, i, j, s, opt;
  pen pena, penb;
  pair PO[];
  u = 1.5mm;
  s = 7.1u;
  pena = pencircle scaled 1.6u;
  penb = pencircle scaled 0.6u;
  linecap := squared;
  for i=1 upto 4:
    for j=1 upto 3:
      draw ((i-1)*s,(j-1)*s)--
	((i-1)*s,j*s)--
	(i*s,j*s)--
	(i*s,(j-1)*s)--cycle withcolor 0.7white;
    endfor;
  endfor;
  z1 = (0,0);
  z2 = (4,3);
  z4 = (1,2);
  z5 = (2,2);
  z6 = (3,2);
  z7 = (1,1);
  z8 = (2,1);
  forever:
    opt := uniformdeviate(1);
    if opt<0.3333:
      PO3 := (1,3);
    elseif opt<0.6667:
      PO3 := (0,3);
    else:
      PO3 := (0,2);
    fi;
    opt := uniformdeviate(1);
    if opt<0.3333:
      PO9 := z4;
    elseif opt<0.6667:
      PO9 := z5;
    else:
      PO9 := z7;
    fi;
    numeric massadois, massatres;
    PO9 = massadois*(z2-PO9)+massatres*(PO3-PO9);
    exitif (massadois>0) and (massatres>0);
  endfor;
  linecap := rounded;
  z3 = PO3;
  for i=1 upto 3:
    draw z[i] scaled s withpen pena;
  endfor;
  for i=4 upto 8:
    draw z[i] scaled s withpen penb;
  endfor;
  label.urt(decimal(round(massatres*24)) & "kg", z3 scaled s);
  label.urt(decimal(round(massadois*24)) & "kg", z2 scaled s);
  label.urt("24kg", z1 scaled s);
  label.urt(btex A etex, z4 scaled s);
  label.urt(btex B etex, z5 scaled s);
  label.urt(btex C etex, z6 scaled s);
  label.urt(btex D etex, z7 scaled s);
  label.urt(btex E etex, z8 scaled s);
  draw fullcircle scaled 1.6u shifted (PO9 scaled s) withcolor red;
endfig;

verbatimtex \end{document} etex

end;
