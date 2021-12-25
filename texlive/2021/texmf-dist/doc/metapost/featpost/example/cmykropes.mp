% cmykropes.mp
% L. Nobre G. 
% 2006

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

verbatimtex
\documentclass{article} 
\usepackage{palatino,mathpazo} 
\begin{document} 
etex

Spread := 39;
f := 3*(3,5,4);

color ki[]t[];
numeric inum, icen, bnum, knum, frame, keytime, maxtetradimstep, atime;
inum = 18;
icen = 9;
bnum = inum*8;
knum = 5;
maxtetradimstep = 0.2;
atime = 4.5;

ki0t0 = (-6.6,2.0,0.3);
ki1t0 = (-5.3,1.5,0.0);
ki2t0 = (-2.3,1.0,-0.4);
ki3t0 = (-0.5,0.6,0.8);
ki4t0 = (0.3,-0.2,-0.4);
ki5t0 = (0.9,-0.2,-1.1);
ki6t0 = (0.9,1.8,-0.6);
ki7t0 = (-0.4,1.8,0.4);
ki8t0 = (-0.4,0.3,-0.2);
ki9t0 = (0.3,-6.0,0.1);
ki10t0 = (5.7,-5.5,0.6);
ki11t0 = (5.8,0.9,0.6);
ki12t0 = (1.6,1.2,0.4);
ki13t0 = (0.1,1.1,-0.6);
ki14t0 = (-2.4,1.4,0.3);
ki15t0 = (-2.4,0.0,0.1);
ki16t0 = (0.2,0.3,-0.9);
ki17t0 = (1.6,0.3,0.3);
ki18t0 = (3.9,0.4,0.4);

ki0t1 = (-6.6,2.0,0.3);
ki1t1 = (-5.3,1.5,0.0);
ki2t1 = (-2.3,1,-0.4);
ki3t1 = (-0.5,0.6,0.8);
ki4t1 = (0.3,-0.2,-0.4);
ki5t1 = (0.9,-0.2,-1.1);
ki6t1 = (0.9,1.8,-0.6);
ki7t1 = (-0.4,1.8,0.4);
ki8t1 = (-0.4,0.3,-0.2);
ki9t1 = (0.3,-6.0,0.1);
ki10t1 = (5.7,-5.5,0.6);
ki11t1 = (5.8,0.9,0.6);
ki12t1 = (1.6,1.2,0.4);
ki13t1 = (0.1,1.1,-0.6);
ki14t1 = (-2.4,1.8,0.3);
ki15t1 = (-5.4,2.0,0.3);
ki16t1 = (-7.2,2.2,0.3);
ki17t1 = (-8.6,2.4,0.3);
ki18t1 = (-9.9,2.8,0.4);

ki0t2 = (-6.6,2,0.3);
ki1t2 = (-5.3,1.5,0.0);
ki2t2 = (-2.3,1.0,-0.4);
ki3t2 = (-0.5,0.6,0.8);
ki4t2 = (0.3,-0.2,-0.4);
ki5t2 = (0.9,-0.2,-1.1);
ki6t2 = (0.9,1.8,-0.6);
ki7t2 = (-0.4,1.8,0.4);
ki8t2 = (-0.4,0.3,-0.2);
ki9t2 = (0.3,-6.0,0.1);
ki10t2 = (5.7,-5.5,0.6);
ki11t2 = (5.8,0.9,0.6);
ki12t2 = (1.6,1.2,0.4);
ki13t2 = (0.1,1.1,-0.6);
ki14t2 = (-2.4,1.8,0.3);
ki15t2 = (-5.4,1.4,0.3);
ki16t2 = (-5.2,1,0.4);
ki17t2 = (-3.6,0.8,0.4);
ki18t2 = (-2.4,0.8,0.4);

ki0t3 = (-6.6,2,0.3);
ki1t3 = (-5.3,1.5,0.0);
ki2t3 = (-2.3,1.0,-0.4);
ki3t3 = (-0.5,0.6,0.8);
ki4t3 = (0.3,-0.2,-0.4);
ki5t3 = (0.9,-0.2,-1.1);
ki6t3 = (0.9,1.8,-0.6);
ki7t3 = (-0.4,1.8,0.4);
ki8t3 = (-0.4,0.3,-0.2);
ki9t3 = (0.3,-6.0,0.1);
ki10t3 = (5.7,-5.5,0.6);
ki11t3 = (5.8,0.9,0.6);
ki12t3 = (1.6,1.2,0.4);
ki13t3 = (0.1,1.1,-0.6);
ki14t3 = (-2.4,1.8,0.3);
ki15t3 = (-3.4,1.4,0.3);
ki16t3 = (-2.4,0.8,0.3);
ki17t3 = (0.1,0.8,-0.6);
ki18t3 = (1.6,0.8,0.4);

ki0t4 = (-6.6,2,0.3);
ki1t4 = (-5.3,1.5,0.0);
ki2t4 = (-2.3,1.0,-0.4);
ki3t4 = (-0.5,0.6,0.8);
ki4t4 = (0.3,-0.2,-0.4);
ki5t4 = (0.9,-0.2,-1.1);
ki6t4 = (0.9,1.8,-0.6);
ki7t4 = (-0.4,1.8,0.4);
ki8t4 = (-0.4,0.3,-0.2);
ki9t4 = (0.3,-6.0,0.1);
ki10t4 = (5.7,-5.5,0.6);
ki11t4 = (5.8,0.9,0.6);
ki12t4 = (1.6,1.2,0.4);
ki13t4 = (0.1,1.1,-0.6);
ki14t4 = (-2.4,1.8,0.3);
ki15t4 = (-2.4,0.8,0.3);
ki16t4 = (0.1,0.8,-0.6);
ki17t4 = (1.6,0.8,0.4);
ki18t4 = (5.0,0.8,0.4);

ki0t5 = (-6.6,2,0.3);
ki1t5 = (-5.3,1.5,0.0);
ki2t5 = (-2.3,1.0,-0.4);
ki3t5 = (-0.5,0.6,0.8);
ki4t5 = (0.3,-0.2,-0.4);
ki5t5 = (0.9,-0.2,-1.1);
ki6t5 = (0.9,1.8,-0.6);
ki7t5 = (-0.4,1.8,0.4);
ki8t5 = (-0.4,0.3,-0.2);
ki9t5 = (0.3,-6.0,0.1);
ki10t5 = (5.7,-5.5,0.6);
ki11t5 = (5.8,0.9,0.6);
ki12t5 = (1.6,1.2,0.4);
ki13t5 = (0.1,1.1,-0.6);
ki14t5 = (-2.4,1.3,0.3);
ki15t5 = (-2.4,0.8,0.3);
ki16t5 = (0.1,0.8,-0.6);
ki17t5 = (1.6,0.8,0.4);
ki18t5 = (5.0,0.8,0.4);


% Hermite polynomials

vardef Pne( expr pr ) = ( (2*pr-3)*(pr**2)+1 ) enddef;
vardef Pwo( expr pr ) = ( (3-2*pr)*(pr**2)   ) enddef;
vardef Tne( expr pr ) = ( pr*(pr*(pr-2)+1)   ) enddef;
vardef Two( expr pr ) = ( (pr-1)*(pr**2)     ) enddef;

% Derivatives

vardef Pone( expr pr ) = ( 6*pr*(pr-1)   ) enddef;
vardef Ptwo( expr pr ) = ( 6*pr*(1-pr)   ) enddef;
vardef Tone( expr pr ) = ( pr*(3*pr-4)+1 ) enddef;
vardef Ttwo( expr pr ) = ( pr*(3*pr-2)   ) enddef;

vardef Wei( expr pre, nex ) =
  ( 0.5*( nex - pre ) )
enddef;

vardef Hermite( expr pr, zer, one, two, thr ) =
  ( Pne(pr)*one+Pwo(pr)*two+Tne(pr)*Wei(zer,two)+Two(pr)*Wei(one,thr) )
enddef;

vardef HermiteRate( expr pr, zer, one, two, thr ) =
  ( Pone(pr)*one+Ptwo(pr)*two+Tone(pr)*Wei(zer,two)+Ttwo(pr)*Wei(one,thr) )
enddef;

vardef RopeThisWay( expr Sens, Istart, Bnum )( text KP ) =
  numeric pstep, p, ap, bray;
  numeric istart, ifinis, thfac, vhfac, sf;
  color actp, newpoi;  
  pstep = 0.002;
  bray = 0.12;
  thfac = 15.5;
  vhfac = 0.8;
  if Sens: sf = 1; else: sf = -1; fi;
  istart = Istart;
  ifinis = Istart+sf;
  ap := 0;
  forever:
    ap := ap+pstep;
    actp := Hermite(ap,KP[istart-sf],KP[istart],KP[ifinis],KP[ifinis+sf]);
    exitif conorm(actp-KP[istart])>bray;
  endfor;
  getready("signalvertex("
      & cstr(actp)
      & "," & decimal(thfac)
      & ",(0,0,0));signalvertex("
      & cstr(actp)
      & "," & decimal(thfac*vhfac)
      & ",(1,1,1));", actp);
  p = 1;
  forever:
    ap := ap+pstep;
    if ap>1:
      ap := ap-1;
      istart := ifinis;
      ifinis := istart+sf;
    fi;
    exitif (ifinis>inum-1) or (p>Bnum) or (ifinis<1);
    newpoi := Hermite(ap,KP[istart-sf],KP[istart],KP[ifinis],KP[ifinis+sf]);
    if conorm(newpoi-actp)>2*bray:
      actp := newpoi;
      getready("signalvertex("
	  & cstr(actp)
	  & "," & decimal(thfac)
	  & ",(0,0,0));signalvertex("
	  & cstr(actp)
	  & "," & decimal(thfac*vhfac)
	  & ",(1,1,1));", actp);
      p := incr( p );
    fi;
  endfor;
enddef;

def RopeBothWays( expr Icen, Tnum )( text VKP ) =
  begingroup
    save bnumpos, bnumneg;
    numeric bnumpos, bnumneg;
    bnumpos = floor(Tnum*(inum-Icen)/inum);
    bnumneg = Tnum-bnumpos;
    RopeThisWay( true, Icen, bnumpos, VKP );
    RopeThisWay( false, Icen, bnumneg, VKP )
  endgroup
enddef;

beginfig(1);
  MalcomX := false;
  numeric i, a, b, c, d, rate, ratenorm, counter;
  cmykcolor ca, cb, cc, cd, aux;
  color kp[];
  frame = 0;
  keytime = 0;
  counter = 0;
  forever:
    a := keytime-1;
    if a<0: a := 0; fi;
    b := keytime;
    c := keytime+1;
    d := keytime+2;
    if d>knum: d := knum; fi;
    ratenorm := 0;
    for i=0 upto inum:
      ca := makecmyk( ki[i]t[a], a );
      cb := makecmyk( ki[i]t[b], b );
      cc := makecmyk( ki[i]t[c], c );
      cd := makecmyk( ki[i]t[d], d );
      aux := HermiteRate( frame, ca, cb, cc, cd );
      rate := cmyknorm( aux );
      if rate>ratenorm: ratenorm := rate; fi;
    endfor;
    frame := frame + W( aux )*maxtetradimstep/ratenorm;
    if frame>1:
      frame := frame-1;
      keytime := keytime+1;
    fi;
    exitif (keytime>knum-1) or (keytime+frame>atime);
  endfor;  
  a := keytime-1;
  if a<0: a := 0; fi;
  b := keytime;
  c := keytime+1;
  d := keytime+2;
  if d>knum: d := knum; fi;
  for i=0 upto inum:
    kp[i] := Hermite( frame, ki[i]t[a], ki[i]t[b], ki[i]t[c], ki[i]t[d] );
  endfor;
%  counter := incr( counter );
%  beginfig(counter);
    setthearena( 15, 14 );
    RopeBothWays( icen, bnum, kp );
    doitnow;
%  endfig;
  endfig;
  
verbatimtex 
\end{document} 
etex

end;
