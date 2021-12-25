% J. Schwaiger with help of
% FEATPOST and the macropackages
% la3Dmacros.mp and hlr3Dmacros.mp
% by L. Nobre G.
% 2003

% Figure (2) retouched in 2004 by L. Nobre G.

input featpost3Dplus2D;

verbatimtex
\documentclass[a4paper,12pt]{article}
\usepackage[T1]{fontenc}
\usepackage[latin1]{inputenc}
\usepackage{euler}
\usepackage{graphicx}
\begin{document}
etex

% for trigonometric functions when angles are measured in rad
%    
    def sinrad( expr t)=sind(t/3.145*180) enddef;
    def cosrad( expr t)=cosd(t/3.145*180) enddef;

    def parsurface( expr tc,sc )=
      (%
	(4 cosrad(tc) - sinrad(sc)*cosrad(tc),%
	  4 sinrad(tc) - sinrad(sc)*sinrad(tc),%
	  cosrad(sc))%
	)% A torus
    enddef;

    def parsurfacetwo( expr tc,sc )=
      (%
	(4+4cosrad(tc) - sinrad(sc)*cosrad(tc),%
	  cosrad(sc),
	  4sinrad(tc) - sinrad(sc)*sinrad(tc)%
	  )%
	)% Another torus
    enddef;

    def parsurfacethree( expr tc,sc )=
% 0<=sc<=2pi,-2<=tc<=2
%xmin=-1.4, xmax=1.4
%ymin=-0.7, ymax=2
%zmin=-2.5,zmax=2.5
      (%
	(1.414213562*cosrad(sc)*cosrad(tc)
	  -2*sinrad(sc)*sinrad(tc))/(2*(abs(tc) + 1)) + sinrad(tc),%
	cosrad(tc) - (1.414213562*cosrad(sc)*sinrad(tc)
	  + 2*sinrad(sc)*cosrad(tc))/(2*(abs(tc) + 1)),% 
	tc - 2*cosrad(sc)/(2*(abs(tc) + 1))
	)%
% A space tube governed by a helix
    enddef;

def zsurface( expr xc, yc ) =
	( xc, yc, cosd(xc*57)*cosd(yc*57)+4*mexp(-(xc**2+yc**2)*6.4) )
enddef;

f:=(8,1.5,1.5);
Spread := 127;


beginfig(1); 
  numeric Ns, Nt, lows, higs, lowt, higt, lowx, highx, lowy, highy,
  facz, lowz, higz, base, topo; 
  path chair;
  % for torus
  f:=(4,1.5,1.5);
  Spread:=150;
  Ns = 16;
  Nt=  16;
  lows := 3.145*0.25;
  higs := 3.145*1.75;
  lowt = -3.145;
  higt = 3.145;
  facz=1;
  lowz = -5;
  higz = 5;
  lowx:=-5;
  lowy:=-5;
  higx:=10;
  higy:=7;
%    
  NF:=0; %Initialize number of faces;
%
  partrimesh( Ns,Nt,lowt,higt,lows,higs,
                    lowx,higx,lowy,higy,lowz,higz,facz,parsurfacetwo );
  lows := 3.145*0.5;
  higs := 2.25*3.145; % new parameters and new surface
  partrimesh( Ns,Nt,lowt,higt,lows,higs,
                    lowx,higx,lowy,higy,lowz,higz,facz,parsurface );
  base = -1/facz;
  topo = 1/facz;
  chair = rp((-1,-1,base))--rp((1,-1,base))--rp((1,1,base))
  --rp((-1,1,base))--rp((-1,-1,base))--rp((-1,-1,topo))
  --rp((-1,1,topo))--rp((-1,1,base))--cycle;
  draw_invisible( false,true,(blue+0.5green),0.6white);
  draw bbox currentpicture;
endfig;
  
beginfig(2); 
  numeric Ns, Nt, lows, higs, lowt, higt, lowx, highx, lowy, highy,
  facz, lowz, higz, base, topo; 
  path chair;

  f := (6,0.5,3); 
  Spread:=150;

  Ns = 19;
  Nt=  19;
  lows := -3.145;
  higs := 3.145;
  lowt = -3.145;
  higt = 3.145;
  facz = 2;
  lowz = -2;
  higz = 2;
  lowx:=-5;
  lowy:=-5;
  higx:=5;
  higy:=5;
  NF:=0; %Initialize number of faces;
%
  partrimesh( Ns,Nt,lowt,higt,lows,higs,
                    lowx,higx,lowy,higy,lowz,higz,facz,parsurface );

  nx := 27;	
  lowy := -7;
  higy := 7;
  lowx := -7;
  higx := 7;
  facz := 2;
  lowz := -1;
  higz := 5;
  partrimesh( nx,nx,lowx,higx,lowy,higy,
	            lowx,higx,lowy,higy,lowz,higz,facz,zsurface);
  base = -1/facz;
  topo = 1/facz;
  chair = rp((-1,-1,base))--rp((1,-1,base))--rp((1,1,base))
  --rp((-1,1,base))--rp((-1,-1,base))--rp((-1,-1,topo))
  --rp((-1,1,topo))--rp((-1,1,base))--cycle;
  fill chair withcolor 0.9 blue+0.6white;
  draw chair withcolor 0.1white ;
  draw_invisible( false,false,(red+0.32blue),0.6white);
endfig;

verbatimtex \end{document} etex
  
end;
