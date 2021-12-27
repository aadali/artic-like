% roadincline.mp
% L. Nobre G.
% 2004
verbatimtex 
\documentclass{article} 
\usepackage{palatino} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex
beginfig(1);
  numeric u, size, halfs, parcel, bordlarg, bordtotal, backdiaglen, backang;
  numeric backdown, labelmarg, backauxil, backauyil, angray;
  path triborder, backtri, angline;
  picture labelalfa, labelv, labelh, labeld;

  u = 1mm;
  size = 25u;
  bordlarg = 3u;
  bordtotal = 4u;
  backdiaglen = 55u;
  backang = 29;
  backdown = 2u;
  labelmarg = -1u;
  angray = 5u;

  parcel = 0.16667*sqrt(3)*size;
  halfs = 0.5size;
  z1 = (-halfs,-parcel);
  z2 = ( halfs,-parcel);
  z3 = ( 0, 2*parcel );
  triborder = z1--z2--z3--cycle;

  backauxil = backdiaglen*cosd(backang)*0.5;
  backauyil = backdiaglen*sind(backang)*0.5;
  z4 = (-backauxil,-backauyil);
  z5 = ( backauxil,-backauyil);
  z6 = ( backauxil, backauyil);
  backtri = z4--z5--z6--cycle;
  
  angline = halfcircle scaled 2angray shifted z4 cutafter (z4--z6);
  z7 = z4+(angray,0);
  z8 = 0.5[z5,z6];
  z9 = z5-(angray,0);
  z0 = z6-angray*(cosd(backang),sind(backang));

  labelalfa = thelabel.urt( btex $\alpha$ etex, z7 );
  labelv = thelabel.lft( btex $v$ etex, z8 );
  labelh = thelabel.top( btex $h$ etex, z9 );
  labeld = thelabel.bot( btex $d$ etex, z0 );

  filldraw backtri;
  draw angline withcolor background;
  draw labelalfa withcolor background;
  draw labelv withcolor background;
  draw labelh withcolor background;
  draw labeld withcolor background;

  draw triborder withcolor background withpen pencircle scaled bordtotal;
  draw triborder withcolor red withpen pencircle scaled bordlarg;

  label.top( btex X\% etex rotated backang, labelmarg*up );
endfig;
verbatimtex 
\end{document} 
etex
end.
