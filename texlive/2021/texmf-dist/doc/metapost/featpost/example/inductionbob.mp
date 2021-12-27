% inductionbob.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

verbatimtex
\documentclass{article} 
\usepackage{beton,concmath,ccfonts} 
\begin{document} 
etex

numeric altura;
altura = 0.35;

f := (4,6,1);
Spread := 200;
LightSource := f+blue;
ShadowOn := false;

def spirala( expr TheVal ) = 
  begingroup
    numeric par[], cred, cgre, cblu, param[];
    par1 = 4.5;
    par2 = 0.3;
    par3 = 0.2*par1;
    param1 = par1*360*TheVal;
    cred = -par2*sind( param1 );
    cgre = ( TheVal - 1.0 )*par3;
    cblu = par2*cosd( param1 );
    ( (cred,cgre,cblu+altura) )
  endgroup 
enddef;

def spiralb( expr TheVal ) = 
  begingroup
    numeric par[], cred, cgre, cblu, param[];
    par1 = 4.0;
    par2 = 0.3;
    par3 = 0.2*par1;
    param1 = par1*360*TheVal;
    cred = par2*sind( param1 );
    cgre = TheVal*par3;
    cblu = -par2*cosd( param1 );
    ( (cred,cgre,cblu+altura) )
  endgroup 
enddef;

def innerloop( expr TheVal ) = 
  begingroup
    numeric par[], cred, cgre, cblu, param[];
    par1 = 8;
    par2 = 0.2;
    param1 = par1+(360-2*par1)*TheVal;
    cred = par2*cosd( param1 );
    cgre = 0;
    cblu = par2*sind( param1 );
    ( (cred,cgre,cblu+altura) )
  endgroup 
enddef;

beginfig(1);
  numeric xouterdev, xinnerdev, scalecir, scalelab;
  color startp, begp, endp, finisp, ap, bp, cp, dp, ep, fp, gp;
  pair startz, begz, endz, finisz, az, bz, cz, dz, ez, fz, gz;
  path innerpath, outerpath, primarcir, secondcir;
  pen mypen;
  scalecir = 3cm;
  scalelab = 4;
  xouterdev = 1.4;
  xinnerdev = 1.0;
  mypen = pencircle scaled 4pt;
  startp = spirala( 0 );
  finisp = spiralb( 1 );
  begp = innerloop( 0 );
  endp = innerloop( 1 );
  ap = ( xouterdev, Y(startp), Z(startp) );
  bp = ( xouterdev, Y(startp), Z(finisp) );
  dp = ( xouterdev, Y(finisp), Z(finisp) );
  cp = 0.5[bp,dp]; 
  ep = ( xinnerdev, Y(begp), Z(begp) ); 
  fp = ( xinnerdev, Y(endp), Z(endp) );
  gp = 0.5[ep,fp];
  startz = rp( startp );
  finisz = rp( finisp );
  begz = rp( begp );
  endz = rp( endp );
  az = rp( ap );
  bz = rp( bp );
  cz = rp( cp );
  dz = rp( dp );
  ez = rp( ep );
  fz = rp( fp );
  gz = rp( gp );
  emptyline( false, 1, black, 0.7background, 700, 0.8, 2, spirala );
  outerpath = startz--az--bz--dz--finisz;
  primarcir = fullcircle scaled scalecir shifted cz;
  emptyline( false, 0.8, black, 0.4*(red+background), 150, 0.8, 2, innerloop );
  innerpath = begz--ez--fz--endz;
  secondcir = fullcircle scaled scalecir shifted gz;
  draw innerpath withpen mypen;
  unfill secondcir;
  draw secondcir withpen mypen;
  label( btex $V_S$ etex scaled scalelab, gz);
  emptyline( true, 1, black, 0.7background, 700, 0.8, 2, spiralb );
  draw outerpath withpen mypen;
  unfill primarcir;
  draw primarcir withpen mypen;
  label( btex $V_P$ etex scaled scalelab, cz );
endfig;

verbatimtex 
\end{document} 
etex

end;
