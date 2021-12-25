% trigonometry.mp
% L. Nobre G.
% 2002

% %input mp-tool;

verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

def PickColor =
  ( uniformdeviate(1) , uniformdeviate(1), 0 )
enddef;

beginfig(1);
  picture sin, cos, tan, ctg, sec, csc;
  color sincol, coscol, tancol, ctgcol, seccol, csccol;
  pair sinl, cosl, tanl, ctgl, secl, cscl;
  path sinp, cosp, tanp, ctgp, secp;
  numeric Ang, U, len, labelscale;
  path diag, fram, quar;
  pair urc, ulc, lrc, tgc, sic, coc, btg;
  pen basepen, thikpen, thinpen;
  Ang = 21;              % Ang < 45
  U = 1.8in;
  labelscale = 2;
  basepen = pencircle scaled 1.2pt;
  thikpen = pencircle scaled 3.4pt;
  thinpen = pencircle scaled 1.2pt;
  sincol = PickColor;
  coscol = PickColor;
  tancol = PickColor;
  ctgcol = PickColor;
  seccol = PickColor;
  csccol = PickColor;
  len = U * cosd( Ang ) / sind( Ang );
  urc = ( len, U );
  ulc = ( 0, U );
  lrc = ( len, 0 );
  tgc = ( U, U * sind( Ang ) / cosd( Ang ) );
  sic = ( U * cosd( Ang ), U * sind( Ang ) );
  coc = ( U * cosd( Ang ), 0 );
  btg = ( U, 0 );
  diag = origin--urc;
  fram = origin--lrc--urc--ulc--cycle;
  quar = halfcircle cutafter ( origin--up );
  quar := quar scaled 2U;
  draw fram withpen basepen;
  draw quar withpen basepen;
  sinl = 0.5[coc,sic];
  sinp = coc--sic;
  cosl = 0.5[origin,coc];
  cosp = origin--coc;
  tanl = 0.5[btg,tgc];
  tanp = btg--tgc;
  ctgl = 0.5[ulc,urc];
  ctgp = ulc--urc;
  secl = 0.5[origin,sic];
  secp = origin--tgc;
  cscl = 0.5[origin,urc];
  draw diag withcolor csccol withpen thikpen;
  draw secp withcolor seccol withpen thinpen;
  draw ctgp withcolor ctgcol withpen thinpen;
  draw tanp withcolor tancol withpen thinpen;
  draw cosp withcolor coscol withpen thinpen;
  draw sinp withcolor sincol withpen thinpen;
  sin = thelabel.lft( btex sin etex scaled labelscale, sinl );
  cos = thelabel.top( btex cos etex scaled labelscale, cosl );
  tan = thelabel.rt( btex tan etex scaled labelscale, tanl );
  ctg = thelabel.bot( btex cotan etex scaled labelscale, ctgl );
  sec = thelabel.top( btex sec etex scaled labelscale, origin );
  csc = thelabel.top( btex cosec etex scaled labelscale, origin );
  draw sin withcolor sincol;
  draw cos withcolor coscol;
  draw tan withcolor tancol;
  draw ctg withcolor ctgcol;
  draw sec rotated Ang shifted secl withcolor seccol;
  draw csc rotated Ang shifted cscl withcolor csccol;
endfig;

verbatimtex 
\end{document} 
etex

end.