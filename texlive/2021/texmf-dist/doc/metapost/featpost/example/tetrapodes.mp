% tetrapodes.mp
% L. Nobre G.
% 2009

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{palatino}
\usepackage{mathpazo}
\begin{document}
etex

f := (10,2,3);
Spread := 50;

def rotatearoundZ( expr P, A ) =
  ((X(P)*cosd(A)-Y(P)*sind(A),X(P)*sind(A)+Y(P)*cosd(A),Z(P)))
enddef;
  
beginfig(1); 
  numeric margia, margib, margic, raycil, raysph, angul;
  color ca, cb;
  path sph;
  margia = 0.10;
  margib = 1.52;
  margic = 0.05;
  raycil = 0.55;
  raysph = 1.00;
  angul = 30;
  ca := (0,-raysph-margia-raycil,0);
  cb := (-raycil-margic,0,raysph+margib);
  rigorousdisc(0,true,ca,raycil,cb-ca);
  ca := rotatearoundZ((-raysph-margia-raycil,0,0),angul);
  cb := rotatearoundZ((0,raycil+margic,-raysph-margib),angul);
  rigorousdisc(0,true,ca,raycil,cb-ca);
  sph = rigorousfearpath(black,raysph);
  unfill sph;
  draw sph;
  ca := rotatearoundZ((raysph+margia+raycil,0,0),angul);
  cb := rotatearoundZ((0,-raycil-margic,-raysph-margib),angul);
  rigorousdisc(0,true,ca,raycil,cb-ca);
  ca := (0,raysph+margia+raycil,0);
  cb := (raycil+margic,0,raysph+margib);
  rigorousdisc(0,true,ca,raycil,cb-ca);
endfig;

beginfig(2); 
    numeric tro;
    tro = 2;
    V1 := (tro,0,0);
    V2 := (0,tro,0);
    V3 := (-tro,0,0);
    V4 := (0,-tro,0);
    V5 := (1,1,1);
    V6 := (-1,1,1);
    V7 := (-1,-1,1);
    V8 := (1,-1,1);
    V9 := (1,1,-1);
    V10 := (-1,1,-1);
    V11 := (-1,-1,-1);
    V12 := (1,-1,-1);
    V13 := (0,0,tro);
    V14 := (0,0,-tro);
    for i=1 upto 14:
      V[i] := rotatearoundZ(V[i],45);
      V[i] := (0.4*X(V[i]),Y(V[i]),1.5*Z(V[i]));
    endfor;
    makeface1(1,5,13,8);
    makeface2(2,6,13,5);
    makeface3(3,7,13,6);
    makeface4(4,8,13,7);
    makeface5(1,5,2,9);
    makeface6(2,6,3,10);
    makeface7(3,7,4,11);
    makeface8(4,8,1,12);
    makeface9(1,9,14,12);
    makeface10(2,10,14,9);
    makeface11(3,11,14,10);
    makeface12(4,12,14,11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := 2*(5,5,1.5);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
endfig;

def spacedistort( expr Pos, Ref, Fac ) =
  ( Pos + cdotprod( Pos, Ref )*Fac*Ref )
enddef;

beginfig(3);
  f := (15,10,5);
  Spread := 50;
  numeric numsidecells, numlimcorner, i, j, k, cx, cy, cz, faca, facb;
  numeric ar, scal;
  color posa, posb, posc, posd, thisoria, thisorib, thisoric,thatori;
  path ball;
  numsidecells = 7;
  numlimcorner = 2*numsidecells+2;
  faca = 0.0;
  facb = 0.0;
  ar = 1.0;
  scal = 0.14;
  thisoria = N(( 1, 1, 1));
  thisorib = N(( 1,-2, 1));
  thisoric = 0.5*sqrt(3)*(1,1,1);
  thatori = ncrossprod(thisoric,thisorib);
  for i=1 upto numsidecells:
    cx := i-0.5*numsidecells-1;
    for j=1 upto numsidecells:
      cy := j-0.5*numsidecells-1;
      for k=1 upto numsidecells:
	cz := k-0.5*numsidecells-1;
	if i+j+k < numlimcorner+1:
	  posa := spacedistort( (cx,cy,cz), thisoria, faca );
	  posa := spacedistort( posa, thisorib, facb );
	  posb := spacedistort( (cx-0.5,cy+0.5,cz), thisoria, faca );
	  posb := spacedistort( posb, thisorib, facb );
	  posc := spacedistort( (cx-0.5,cy,cz+0.5), thisoria, faca );
	  posc := spacedistort( posc, thisorib, facb );
	  posd := spacedistort( (cx,cy+0.5,cz+0.5), thisoria, faca );
	  posd := spacedistort( posd, thisorib, facb );
% 	  whatisthis(posa,scal*ar*thatori,scal*thisorib,1,1);
% 	  whatisthis(posb,scal*ar*thatori,scal*thisorib,1,1);
% 	  whatisthis(posc,scal*ar*thatori,scal*thisorib,1,1);
% 	  whatisthis(posd,scal*ar*thatori,scal*thisorib,1,1);
	  ball := rigorousfearpath( posa, scal ); unfill ball; draw ball;
	  ball := rigorousfearpath( posb, scal ); unfill ball; draw ball;
	  ball := rigorousfearpath( posc, scal ); unfill ball; draw ball;
	  ball := rigorousfearpath( posd, scal ); unfill ball; draw ball;
	fi;
      endfor;
    endfor;
  endfor;
  
endfig;

end;
 
