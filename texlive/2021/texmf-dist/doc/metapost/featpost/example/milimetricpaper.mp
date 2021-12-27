% milimetricpaper.mp
% L. Nobre G.
% 2014

prologues := 1;

beginfig(1);
  numeric u, xsize, ysize, i, ticsize;
  numeric minstep, majstep, bigstep;
  pen minpen, majpen, bigpen;
  color mincolor, majcolor, bigcolor;
  u = 1mm;
  xsize = 15cm;
  ysize = 11cm;
  ticsize = 3u;
  minstep = 2u;
  majstep = 10u;
  bigstep = 20u;
  minpen = pencircle scaled 0pt;
  majpen = pencircle scaled 0.4pt;
  bigpen = pencircle scaled 0.99pt;
  mincolor = 0.65white;
  majcolor = 0.35white;
  bigcolor = 0.15white;
  for i=0 step minstep until xsize:
    draw (i,0)--(i,ysize) withpen minpen withcolor mincolor;
  endfor;
  for i=0 step minstep until ysize:
    draw (0,i)--(xsize,i) withpen minpen withcolor mincolor;
  endfor;
  for i=0 step bigstep until xsize:
    draw (i,-ticsize)--(i,ysize) withpen bigpen withcolor bigcolor;
  endfor;
  for i=0 step bigstep until ysize:
    draw (-ticsize,i)--(xsize,i) withpen bigpen withcolor bigcolor;
  endfor;
  for i=0 step majstep until xsize:
    draw (i,0)--(i,ysize) withpen majpen withcolor majcolor;
  endfor;
  for i=0 step majstep until ysize:
    draw (0,i)--(xsize,i) withpen majpen withcolor majcolor;
  endfor;
endfig;

beginfig(2);
  drawarrow origin--(xsize,0) withpen bigpen;
  for i=majstep step majstep until xsize-2majstep:
    draw (i,-minstep)--(i,minstep) withpen majpen;
  endfor;
endfig;

end.
