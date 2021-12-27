
prologues:=1;

numeric u;
1u = 1mm;


beginfig (1);

 a = 40u;
 g = 90;

 z1 = right rotated g scaled a;
 z2 = z1 rotated 120;
 z3 = z1 rotated 240;

 path c[];
 c1 = fullcircle rotated g scaled 1/2a;
 for k = 2 upto 4:
   c[k] = c1 scaled k;
 endfor

 p = 2q;
 p + q = 1/4a;
 d = 2q;

 path r[], s[], s[]';
 for k = 1 upto 3:
   r[k] = (0,0)--z[k];
   s[k] = r[k] shifted (down scaled (p/2+d) rotated angle z[k]);
   s[k]' = r[k] shifted (up scaled (p/2+d) rotated angle z[k]);
 endfor

 pen rpen, cpen;
 cpen := pencircle scaled p;
 rpen := pencircle scaled p;

 pickup cpen;
 filldraw c1 withcolor 0.75 red;

 linecap := rounded;
 draw c2 cutbefore s1' cutafter s1;
 draw c3 cutbefore s1' cutafter s2;
 draw c3 cutbefore s2' cutafter s1;
 draw c4 cutbefore s1' cutafter s2;
 draw c4 cutbefore s2' cutafter s3;
 draw c4 cutbefore s3' cutafter s1;

endfig;


end