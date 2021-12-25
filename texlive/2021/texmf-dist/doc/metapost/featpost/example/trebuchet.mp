% trebuchet.mp
% L. Nobre G.
% 2014

%prologues := 1;
%outputformat := "png";
%outputformatoptions := "format=rgb antialias=best";
%hppp := 0.25;
%vppp := 0.25;

beginfig(1);
  numeric i, s, a, u, h, t;
  pair pi, pj, pk;
  u = 9mm;
  h = 5u;
  a = 0.039u;
  s = 5;
  draw (0,-h) for i=s-90 step s until 90:
      ..(h*(2*cosd(i),sind(i)))
  endfor --cycle;
  draw origin--(2h,0);
  t = 0;
  forever:
    t := incr(t);
    i := h-0.5*a*t**2;
    exitif i<-h;
    pi := (0,i);
    pj := (h+-+i,0);
    pk := 2[pi,pj];
    if i<0:
      draw pj withpen pencircle scaled 0.17u;
    fi;
    draw pi--pk withcolor 0.5white;
  endfor;
endfig;
end.
