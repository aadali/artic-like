% simplelens.mp
% L. Nobre G.
% 2013

prologues := 1;

numeric specialdist, otherdist;
specialdist = 2.5cm;
color TableC[];
TableC0 := (0.65,0.61,0.49);    % grey        %% G N U P L O T
TableC1 := 0.77red;             % red         %% 
TableC2 := ( 0.2, 0.2, 1.0 );   % blue        %%    colors 
TableC3 := ( 1.0, 0.7, 0.0 );   % orange      %%
TableC4 := 0.85green;           % pale green  %%
TableC5 := 0.90*(red+blue);     % magenta     %%
TableC6 := 0.85*(green+blue);   % cyan        %%
TableC7 := 0.85*(red+green);    % yellow      %%

vardef lensrefractray(expr s,a,p,q,rir,c)=
  save tI,tn,ia,ib,I,J,K,sib;
  pair I,J,K;
  draw s--(s-2cm*dir(a)) withcolor TableC[c];
  tI = xpart(p intersectiontimes (s--(s+15cm*dir(a))));
  if tI>=0:
    I=point tI of p;
    draw s--I withcolor TableC[c];
    tn=angle(direction tI of p)+90;
    ia=tn-angle(s-I);
    sib=sind(ia)/rir;
    ib=-angle(1+-+sib,sib);    
    J=I+dir(tn+180+ib)*2cm;
    tI := xpart(q intersectiontimes (I--J));
    if tI>=0:
      K:=point tI of q;
      draw I--K withcolor TableC0;
      tn:=angle(direction tI of q)+90;
      ia:=tn-angle(I-K);
      if sind(ia)<1/rir:
	sib:=rir*sind(ia);
	ib:=angle(1+-+sib,sib);    
	J:=K-dir(tn-ib)*specialdist;
	draw K--J withcolor TableC[c];
      fi;
    fi;
  fi;
enddef;

beginfig(1);
  numeric rir, ray, wid, dist, ha, dx, dy, bord, cx, cy, cent, a, sa, nr, m;
  numeric refang, colorcounter;
  path vseg, cseg, lens, cell;
  pen fordot;
  color lenscolor, cellcolor;
  refang = 41;
  m = 1.2mm;
  nr = 9;
  rir = 1.2;
  ray = 2.3cm;
  wid = 3.9cm;
  dist = 1.1cm;
  bord = 1mm;
  fordot = pencircle scaled 1mm;
  lenscolor = 0.6white;
  cellcolor = 0.4white;
  dy = 0.5*wid;
  dx = ray+-+dy;
  cx = 0.5*bord;
  cy = dy-bord;
  cent = dx-dist;
  ha = angle(dist,cy);
  sa = 2ha/nr;
  z1 = (dx,ray);
  z2 = (dx,-ray);
  z3 = (cent,0);
  z4 = (dx,cy);
  z5 = (dx,-cy);
  z6 = (0,cy);
  z7 = (0,-cy);
  z8 = (ray,cy);
  z9 = (ray,-cy);
  vseg = z2--z1;
  %cseg = halfcircle rotated -90 scaled 2ray;
  cseg = (dx,-dy){dir(refang)}..(ray,0)..{dir(180-refang)}(dx,dy);
  lens = buildcycle(reverse vseg,z7--z9,cseg,z8--z6);
  fill lens withcolor lenscolor;
  cell = ((cx,cy)--(-cx,cy)--(-cx,-cy)--(cx,-cy)--cycle) shifted z3;
  fill cell withcolor cellcolor;
  %draw z3 withpen fordot;
  colorcounter = 1;
  for a=0.5*sa step sa until ha-sa:
    if not (colorcounter<>7):
      colorcounter:= 1;
    else:
      colorcounter := incr(colorcounter);
    fi;
    lensrefractray(z3,a,vseg,cseg,1.9,colorcounter);
    lensrefractray(z3+up*m,a,vseg,cseg,1.9,colorcounter);
    lensrefractray(z3+down*m,a,vseg,cseg,1.9,colorcounter);
    if not (colorcounter<>7):
      colorcounter:= 0;
    else:
      colorcounter := incr(colorcounter);
    fi;
    lensrefractray(z3,-a,vseg,cseg,1.9,colorcounter);
    lensrefractray(z3+up*m,-a,vseg,cseg,1.9,colorcounter);
    lensrefractray(z3+down*m,-a,vseg,cseg,1.9,colorcounter);
  endfor;
  otherdist = ray+specialdist-0.5cm;
  clip currentpicture to (-otherdist,-cy)--(-otherdist,cy)--(otherdist,cy)--(otherdist,-cy)--cycle; 
endfig;

end.
  
  