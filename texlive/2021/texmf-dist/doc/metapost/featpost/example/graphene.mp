% graphene.mp
% L. Nobre G.
% 2011

prologues := 1;
  
input featpost3Dplus2D;

f := (9,13,7);
Spread := 15;

beginfig(1);
  numeric s, c, rs, rc, i, j;
  boolean isodd;
  color ua, ub, uc, ta, tb;
  path t;
  s = 0.5;
  c = sqrt(3)/2;
  rs = 0.23;
  rc = 0.15;
  ua = (-s,c,0);
  ub = (s,c,0);
  uc = (1,0,0);
  fill VGAborder withcolor 0.65white;
  isodd = true;
  j = 1;
  ta = (-4.5,-5,0);
  t = rigorousfearpath( ta, rs );
  unfill t; draw t;
  forever:
    if isodd:
      tb := ta + ua;
      rigorousdisc( 0, true, ta+rs*ua, rc, ua*(1-2*rs));
    else:
      tb := ta + ub;
      rigorousdisc( 0, true, ta+rs*ub, rc, ub*(1-2*rs));
    fi;
    t := rigorousfearpath( tb, rs );
    unfill t; draw t;
    isodd := not isodd;
    ta := tb;
    j := incr( j );
    exitif j>11;
  endfor;
  for i=1 upto 6:
    for j=1 upto 6: 
      if isodd:
	ta := ((i-1)*1.5-4.5,(j-1)*2*c-5+c,0);
      else:
	ta := ((i-1)*1.5-4.5,(j-1)*2*c-5,0);
      fi;
      tb := ta + uc;
      rigorousdisc( 0, true, ta+rs*uc, rc, uc*(1-2*rs));
    endfor;
    j := 1;
    if isodd:
      ta := ((i-1)*1.5-3,j-6,0);
    else:
      ta := ((i-1)*1.5-3.5,j-6,0);
    fi;    
    t := rigorousfearpath( ta, rs );
    unfill t; draw t;
    forever:
      if isodd:
	tb := ta + ua;
	rigorousdisc( 0, true, ta+rs*ua, rc, ua*(1-2*rs));
      else:
	tb := ta + ub;
	rigorousdisc( 0, true, ta+rs*ub, rc, ub*(1-2*rs));
      fi;
      t := rigorousfearpath( tb, rs );
      unfill t; draw t;
      isodd := not isodd;
      ta := tb;
      j := incr( j );
      exitif j>11;
    endfor;
  endfor;
endfig;

end.