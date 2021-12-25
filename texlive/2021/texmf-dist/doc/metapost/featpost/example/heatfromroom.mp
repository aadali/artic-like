% heatfromroom.mp
% L. Nobre G. 
% 2006

verbatimtex 
\documentclass{article} 
\usepackage{palatino,mathpazo} 
\begin{document} 
etex

beginfig(1);
  numeric l, gang;
  path longp, shorp, roomp;
  picture sp[];
  pen isol;
  color mygrey;
  mygrey = 0.75white;
  gang = -angle( (3.5,1) );
  l = 1cm;
  isol = pensquare scaled 2.2pt rotated gang;
  z1  = (0,-l) rotated gang;
  z2  = (l,-l) rotated gang;
  z3  = (l,2*l) rotated gang;
  z4  = (0,2*l) rotated gang;
  z5  = (-3*l,l) rotated gang;
  z6  = (-3*l,0) rotated gang;
  z7  = (3*l,0) rotated gang;
  z8  = (3*l,l) rotated gang;
  z9  = (l,l) rotated gang;
  z10 = (0,l) rotated gang;
  z11 = (l,0) rotated gang;
  z12 = (2*l,0.5*l) rotated gang;
  z13 = (0.5*l,1.5*l) rotated gang;
  z14 = (-1.5*l,0.5*l) rotated gang;
  z15 = (0.5*l,-0.5*l) rotated gang;
  z16 = (l,-0.5*l) rotated gang;
  z17 = (2*l,l) rotated gang;
  z18 = (0,1.5*l) rotated gang;
  z19 = (-1.5*l,0) rotated gang;
  z20 = (0.5*l,-l) rotated gang;
  z21 = (3*l,0.5*l) rotated gang;
  z22 = (0.5*l,2*l) rotated gang;
  z23 = (-3*l,0.5*l) rotated gang;
  z24 = (0.5*l,0.5*l) rotated gang;
  longp = z5--z6--z7--z8--cycle;
  shorp = z1--z2--z3--z4--cycle;
  roomp = origin--z11--z9--z10--cycle;
  fill longp withcolor mygrey;
  fill shorp withcolor mygrey;
  unfill roomp;
  draw longp;
  draw shorp;
  drawoptions( withpen isol );
  draw z1--origin--z6;
  draw z5--z10--z4;
  draw z3--z9--z8;
  draw z7--z11--z2;
  drawoptions();
  sp12 = btex C etex;
  label( sp12, z12 );
  sp13 = btex B etex;
  label( sp13, z13 );
  sp14 = btex A etex;
  label( sp14, z14 );
  sp15 = btex D etex;
  label( sp15, z15 );
  sp16 = btex $L$ etex;
  label.rt( sp16, z16 );
  sp17 = btex $2L$ etex;
  label.top( sp17, z17 );
  sp18 = btex $L$ etex;
  label.lft( sp18, z18 );
  sp19 = btex $3L$ etex;
  label.bot( sp19, z19 );
  sp20 = btex $5^\circ$C etex;
  label.bot( sp20, z20 );
  sp21 = btex -$25^\circ$C etex;
  label.rt( sp21, z21 );
  sp22 = btex -$5^\circ$C etex;
  label.top( sp22, z22 );
  sp23 = btex -$10^\circ$C etex;
  label.lft( sp23, z23 );
  sp24 = btex $5^\circ$C etex;
  label( sp24, z24 );
endfig;

verbatimtex \end{document} etex

end;

