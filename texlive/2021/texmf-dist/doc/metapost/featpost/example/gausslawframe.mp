% gausslawframe.mp 
% adapted from an ANEEMATE example
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

numeric ParameterA, ParameterB, ParameterC, ParameterD, ParameterE;
ParameterA = 0; % 
ParameterB = 0.02;  % escala geral
ParameterC = 0.8; % avanco da superficie
ParameterD = 1; % deslocamento do viewcentr
ParameterE = 3.1; % capitulo

verbatimtex
\documentclass{article} 
\usepackage{mathpazo}
\newfont{\featfont}{logod10 scaled 500}
\begin{document} 
etex

def sfsur( expr H, V ) =
  ( (cosd(H)*cosd(V),sind(H)*cosd(V),sind(V)) )
enddef;

primarydef a joi b = a & "," & b enddef;
primarydef a empstr b = a & ditto & ditto & b enddef;

beginfig(1);
  color auxb, auxc, pa, pb, pc, pd;
  numeric aux, i, j, cou, lim, max, numc, numofdiarr, radratio;
  numeric astepv, asteph, radiub, arch, starter, coverlim, fieldlen;
  numeric refav, refah, actad, minad, actav, actah, cop, lscale, maxlscale;
  pen thickp;
  string a, b, auxs;
  picture AM, DPHI;
  path p[];
  pair auxp, lasp, lpoi, lvec;;
  f 		:= (2,3.5,4.6);
  Spread 	:= ParameterB[9,880];
  fieldlen      = 0.95;
  radius	= 1;
  radratio      = 3.5;
  numc		= 9;
  numofdiarr    = 120;
  refah         = 50;
  refav         = 30;
  actad         = ((refav+90)/180) ++ (refah/180);
  minad         = actad;
  actah         = refah;
  actav         = refav;
  cop           = 0;
  thickp	= pencircle scaled (ParameterB[0.2mm,1.2mm]);
  maxlscale     = 3;
  lscale	= ParameterB[1,maxlscale];
  max           = ceiling( radratio*numc*0.5 );
  astepv        = 180/max;
  radiub        = radius*radratio;
  arch          = radiub*astepv;
  coverlim      = ParameterC[3*astepv-90,89];
  viewcentr     := ParameterD[black,0.5*(radiub+fieldlen)*sfsur(actah,actav)];
  pickup thickp;
  getready( "tropicalglobe( numc, black, radius, blue );", black );
  TDAtiplen := 0.15;    
  TDAhalftipbase := 0.06;
  TDAhalfthick := 0.03;
  cou = 0;
  for j=1.5*astepv-90 step astepv until coverlim:
    asteph := arch/(radiub*cosd(j-0.5*astepv));
    lim := ceiling(360/asteph);
    asteph := 360/lim;
    starter := uniformdeviate( asteph );
    for i=starter-180 step asteph until 180:
      cou := incr( cou );
      pa := radiub*sfsur(i-asteph,j-astepv);
      pb := radiub*sfsur(i,j-astepv);
      pc := radiub*sfsur(i,j);
      pd := radiub*sfsur(i-asteph,j);
      actad := ((i-0.5*asteph-refah)/180) ++ ((j-0.5*astepv-refav)/180);
      if actad < minad:
	minad := actad;
	actav := j-0.5*astepv;
	actah := i-0.5*asteph;
	cop := cou;
      fi;
      auxc := radiub*sfsur( i-0.5*asteph, j-0.5*astepv );
      p[cou] = rp(pa)--rp(pb)--rp(pc)--rp(pd)--cycle;
      auxs := decimal( cou );
      getready( "unfill p" & auxs & "; draw p" & auxs & ";", auxc ); 
      a := cstr( auxc );
      b := cstr( auxc+fieldlen*N(auxc) );
      getready( "tdarrow(" & a joi b & ");", 1.1*auxc );
    endfor;
  endfor;
  doitnow;
  draw VGAborder;
  AM = btex {\featfont ANEEMATE} etex;
  label.ulft( AM rotated 90, point 1 of VGAborder);
  unfill p[cop];
  draw p[cop];
  auxb := radiub*sfsur(actah,actav);
  auxc := (radiub+fieldlen)*sfsur(actah,actav);
  tdarrow( auxb, auxc );
  drawoptions( withcolor blue );
  label.lft( btex $\vec{E}$ etex scaled lscale, rp( auxc ) );
  drawoptions( withpen thickp withcolor red );
  auxc := (radiub+0.5*fieldlen)*sfsur(actah,actav);
  drawarrow rp(auxb)--rp(auxc);
  label.rt( btex $d\vec{A}$ etex scaled lscale, rp( auxc ) );
  drawoptions( withcolor blue );
  if (ParameterE>0.5) and (ParameterE<1.5):
    DPHI = btex $d\Phi=\vec{E}\cdot d\vec{A} = E\,dA$ etex scaled maxlscale;
    label.lrt( DPHI, point 3 of VGAborder);
  fi;
  if (ParameterE>1.5) and (ParameterE<3.5):
    DPHI = btex $\Phi=\int d\Phi
                     =\int \vec{E}\cdot d\vec{A}
		     =\int E\,dA
		     =E\int dA
		     =4\pi R^2 E
		     =\frac{Q}{\mathcal{E}_0}
		     \Leftrightarrow
		     E= \frac{Q}{4\pi R^2\mathcal{E}_0}$ etex scaled maxlscale;
    auxp = ( lrcorner DPHI ) - ( llcorner DPHI ); 
    lasp = ( point 3 of VGAborder ) + auxp;
    lvec = lasp - ( point 2 of VGAborder );
    lpoi = ( point 3 of VGAborder ) - (ParameterE-2)*lvec;
    label.lrt( DPHI, lpoi );
  fi;
  
  produce_vga_border;
endfig;

verbatimtex 
\end{document} 
etex

end;







