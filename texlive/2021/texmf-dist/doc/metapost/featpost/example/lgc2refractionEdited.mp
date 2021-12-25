% lgc2refractionEdited.mp
% lgc means LaTeX Graphics Companion
% 2014

prologues := 1;

vardef refractray(expr s,a,p,rir)=
  save tI,tn,ia,ib,I,J,sib;
  pair I,J;
  tI = xpart(p intersectiontimes (s--(s+30cm*dir(a))));
  if tI>=0:
    I=point tI of p;
    draw s--I;
    tn=angle(direction tI of p)+90;
    ia=tn-angle(s-I);
    sib=sind(ia)/rir;
    ib=-angle(1+-+sib,sib);
    J=I+dir(tn+180+ib)*2cm;
    draw I--J withcolor blue;
  fi;
enddef;

beginfig(1);
  numeric u;
  u=3mm;
  path p;
  p=((5u,-6u)..(4u,-2u)..(4.5u,0)..(5u,5u)) rotated 90;
  for a=40 step 5 until 140: refractray(origin,a,p,1.5); endfor;
  draw p withcolor red;
endfig;

vardef lensrefractray(expr s,a,p,q,rir)=
  save tI,tn,ia,ib,I,J,K,sib;
  pair I,J,K;
  tI = xpart(p intersectiontimes (s--(s+30cm*dir(a))));
  if tI>=0:
    I=point tI of p;
    draw s--I;
    tn=angle(direction tI of p)+90;
    ia=tn-angle(s-I);
    sib=sind(ia)/rir;
    ib=-angle(1+-+sib,sib);    
    J=I+dir(tn+180+ib)*2cm;
    tI := xpart(q intersectiontimes (I--J));
    if tI>=0:
      K:=point tI of q;
      draw I--K withcolor blue;
      tn:=angle(direction tI of q)+90;
      ia:=tn-angle(I-K);
      sib:=rir*sind(ia);
      ib:=angle(1+-+sib,sib);    
      J:=K+dir(tn+ib)*2cm;
      draw K--J;
    fi;
  fi;
enddef;

beginfig(2);
  path q;
  q = reverse (p yscaled (-1)) shifted (up*10.5u);
  fill (p--q--cycle) withcolor red;
  for a=46 step 3 until 120:
    lensrefractray(origin,a,p,q,1.4);
  endfor;
endfig;

vardef reflectray(expr s,a,p,l)=
  save tI,tn,ia,I,J;
  pair I,J;
  tI=xpart(p intersectiontimes (s--(s+30cm*dir(a))));
  if tI>=0:
    I=point tI of p;
    draw s--I;
    tn=angle(direction tI of p)+90;
    ia=tn-angle(s-I);
    J=I+dir(tn+ia)*l;
    draw I--J withcolor blue;
  fi;
enddef;

beginfig(3);
  for a=40 step 5 until 140: reflectray(origin,a,p,2cm); endfor;
  draw p withcolor red;
endfig;

vardef doublereflection(expr s,a,p,q,l,m)=
  save tI,tn,ia,I,J,K;
  pair I,J,K;
  tI=xpart(p intersectiontimes (s--(s+30cm*dir(a))));
  if tI>=0:
    I=point tI of p;
    draw s--I;
    tn=angle(direction tI of p)+90;
    ia=tn-angle(s-I);
    J=I+dir(tn+ia)*l;
    tI := xpart(q intersectiontimes (I--J));
    if tI>=0:
      K:=point tI of q;
      draw I--K withcolor blue;
      tn:=angle(direction tI of q)+90;
      ia:=tn-angle(I-K);
      J:=K+dir(tn+ia)*m;
      draw K--J;
    fi;
  fi;
enddef;

beginfig(4);
  numeric u, ad, len;
  u=3mm;
  ad=83;
  len = 5cm;
  path p,q;
  pair sourc;
  sourc = (-2u,0);
  p=(-5u,0){up}..{right}(3u,9u);
  q=(6u,0)..{dir(ad)}(10u,11u);
  for a=60 step 3 until 120:
    %reflectray(sourc,a,p,len);
    doublereflection(sourc,a,p,q,len,0.3*len);
  endfor;
  draw p withcolor red;
  draw q withcolor red;
endfig;

end.