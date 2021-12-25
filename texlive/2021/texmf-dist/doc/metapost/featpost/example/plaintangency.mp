% plaintangency.mp
% L. Nobre G.
% IYP (2005)

%input mp-tool;

def paircrossprod(expr A, B) = 
  ( (xpart A)*(ypart B) - (xpart B)*(ypart A) )
enddef;

def firsttangencypoint( expr Path, Point, ResolvN ) =
  begingroup
    save auxp, i, cutp, va, vb;
    path auxp;
    numeric i;
    pair cutp, va, vb;
    auxp =
    hide( va := unitvector( point 0 of Path - Point );
      vb := unitvector( direction 0 of Path ); )
    ( paircrossprod( va, vb ), 0 )
    for i=1/ResolvN step 1/ResolvN until length Path:
      hide( va := unitvector( point i of Path - Point );
	vb := unitvector( direction i of Path ); )
      ...( paircrossprod( va, vb ), i )
    endfor;
    cutp = auxp intersectionpoint ( origin--( 0, length Path ) );
    ( point ( ypart cutp ) of Path )
  endgroup
enddef;

beginfig(0);
  numeric u, i;
  u = 5mm;
  pen a, b, c;
  a = pencircle scaled 3pt;
  b = pencircle scaled 5pt;
  c = pencircle scaled 1pt;
  z1 = (1u,1u);
  z2 = (4u,4u);
  z3 = (4u,5u);
  z4 = (3u,5u);
  z5 = (3u,6u);
  z6 = (4u,7u);
  z7 = (6u,1u);
  path cp;
  cp = z1{up}..z2..z3..z4..z5..{up}z6;
  draw cp withpen c;
  for i=1 upto 6:
    draw z[i] withpen a withcolor 0.5*(red+green);
  endfor;
  z8 = firsttangencypoint( cp, z7, 5 );
  draw z7 withpen b withcolor green;
  draw z7--z8 withpen c withcolor blue;
endfig;

beginfig(1);
  numeric u;
  u = 5mm;
  pen a, b, c;
  path auxp;
  numeric i, auxn, yfact, sinfact, res;
  pair cutp, vA, vB;
  yfact = 30;
  sinfact = 40;
  res = 15;
  a = pencircle scaled 3pt;
  b = pencircle scaled 5pt;
  c = pencircle scaled 1pt;
  z1 = (1u,1u);
  z2 = (4u,4u);
  z3 = (4u,5u);
  z4 = (3u,5u);
  z5 = (3u,6u);
  z6 = (4u,7u);
  z7 = (6u,2.8u);
  path cp;
  cp = z1{up}..z2..z3..z4..z5..{up}z6;
  draw z7 withpen b withcolor green;
  draw cp withpen a;
  auxp = hide( vA := unitvector(point 0 of cp - z7); 
               vB := unitvector(direction 0 of cp); )
    ( sinfact*((xpart vA)*(ypart vB) - (xpart vB)*(ypart vA)), 0 )
    for i=1/res step 1/res until length cp:
      hide( vA := unitvector(point i of cp - z7); 
            vB := unitvector(direction i of cp); )
    ...(sinfact*((xpart vA)*(ypart vB)-(xpart vB)*(ypart vA)),i*yfact)
    endfor;
  draw auxp withcolor blue+green;
  draw origin--(sinfact,0) withcolor red;
  draw origin--( 0, yfact*length cp ) withcolor red+green;
  cutp = auxp intersectionpoint ( origin--( 0, yfact*length cp ) );
  draw cutp withpen a;
  auxn = ( ypart cutp )/yfact;
  show auxn;
  z8 = point auxn of cp;
  draw z8 withpen b withcolor green;
  draw z7--1.8[z7,z8] withpen c withcolor blue;
endfig;

end.






