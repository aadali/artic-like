% twoplustwo.mp
% L. Nobre G. 
% 2006

def screenpath( expr CentrePoi, Width, Height, CorneRay ) =
  begingroup
    pair sp[];
    sp0 = CentrePoi - 0.5*( Width*right + Height*up ) + CorneRay*down;
    sp1 = CentrePoi + 0.5*( Width*right - Height*up ) + CorneRay*down;
    sp2 = CentrePoi + 0.5*( Width*right - Height*up ) + CorneRay*right;
    sp3 = CentrePoi + 0.5*( Width*right + Height*up ) + CorneRay*right;
    sp4 = CentrePoi + 0.5*( Width*right + Height*up ) + CorneRay*up;
    sp5 = CentrePoi - 0.5*( Width*right - Height*up ) + CorneRay*up;
    sp6 = CentrePoi - 0.5*( Width*right - Height*up ) + CorneRay*left;
    sp7 = CentrePoi - 0.5*( Width*right + Height*up ) + CorneRay*left;
    ( sp0---sp1..sp2---sp3..sp4---sp5..sp6---sp7..cycle )
  endgroup
enddef;
    
beginfig(1);
  numeric u, i, j, k, l;
  pen stdpen;
  color bgcol, fgcol;
  bgcol = 0.7white;
  fgcol = 0.85*(blue+green);
  stdpen = pencircle scaled 9pt;
  u = 5cm;
  for i=1 upto 5:
    j := 90-72+i*72;
    k := i+5;
    l := k+5;
    z[i] = u*dir(j);
    z[k] = u*sind(18)*dir(j);
    z[l] = u*sind(18)*dir(j+36)/cosd(36);
  endfor;
  path inip, inis, init, iniq, scrbord;
  scrbord = screenpath( origin, 1.7u, 1.1u, 0.7u );
  inip = z5---z2--z12{unitvector(z8-z2)}..{unitvector(z5-z9)}z14--cycle;
%  inis = ( inip ) rotated 72;
%  init = ( ( ( inip ) shifted -z5 ) yscaled -1 ) shifted z4;
%  iniq = ( ( ( ( inip ) shifted -z2 ) yscaled -1 ) rotated 72 ) shifted z4;
%  draw inip withpen pencircle scaled 3pt;
%  draw inis;
%  draw init;
%  draw iniq;
  fill inip withcolor fgcol;
  draw inip withpen stdpen withcolor bgcol;
  picture unitpic;
  unitpic = currentpicture;
  currentpicture := nullpicture;
  fill scrbord withcolor bgcol;
  draw ( ( ( ( unitpic ) shifted -z2 ) yscaled -1 ) rotated 72 ) shifted z4;
  draw unitpic;
  draw ( ( ( unitpic ) shifted -z5 ) yscaled -1 ) shifted z4;
  draw ( unitpic ) rotated 72;
endfig;

end;

