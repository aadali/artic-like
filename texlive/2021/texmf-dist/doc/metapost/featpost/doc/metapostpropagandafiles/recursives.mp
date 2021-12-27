% recursives.mp
% L. Nobre G.
% 2011

prologues := 1;

% coded by Martin Geisler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def koch(expr a, b, n) =
  begingroup
    if n = 0:
      draw a--b;
    else:
      koch(a, 1/3[a,b], n-1);
      koch(1/3[a,b], 2/3[a,b] rotatedabout(1/3[a,b], 60), n-1);
      koch(2/3[a,b] rotatedabout(1/3[a,b], 60), 2/3[a,b], n-1);
      koch(2/3[a,b], b, n-1);
    fi
  endgroup
enddef;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def kochl(expr a, b, n) =
  begingroup
    if n = 0:
      draw a--b;
    else:
      kochl(a, 3/10[a,b], n-1);
      kochl(3/10[a,b], 3/5[a,b] rotatedabout(3/10[a,b], 48.19), n-1);
      kochl(3/5[a,b] rotatedabout(3/10[a,b], 48.19), 7/10[a,b], n-1);
      kochl(7/10[a,b], b, n-1);
    fi
  endgroup
enddef;

beginfig(1);
  numeric u, N;
  u = 10cm;
  N = 5;
  pair a, b, c;
  a = origin;
  b = 0.5*u*(1,sqrt(3));
  c = (u,0);
  koch( a, b, N );
  koch( b, c, N );
  koch( c, a, N );
endfig;

beginfig(2);
  numeric u, N;
  u = 10cm;
  N = 4;
  pair a, b, c;
  a = origin;
  b = 0.5*u*(1,sqrt(3));
  c = (u,0);
  kochl( a, c, N );
  kochl( c, b, N );
  kochl( b, a, N );
endfig;

% from Syracuse
beginfig(3);
  z4 = (u,0);
  z1 = 1/3[origin,z4];
  z3 = 2/3[origin,z4] ;
  z2 - z1 = (z3 - z1) rotated 60 ;
  z4 - origin = (z5 - origin) rotated 60 ;
  def next(expr p) =
    (p scaled 1/3)--(p scaled 1/3 rotated 60 shifted z1)
      --(p scaled 1/3 rotated -60 shifted z2)--(p scaled 1/3 shifted z3)
  enddef;
  path et;
  et = origin{dir(-20)}..{dir(20)}z4;
  for i = 1 upto 3:
    et := next(et) ;
  endfor
  draw et--(et rotated -120 shifted z4)--(et rotated 120 shifted z5);
endfig ;

end.
