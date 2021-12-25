% fis3prex7.mp
% L. Nobre G.
% 2012

input featpost3Dplus2D;

verbatimtex \documentclass{article}\usepackage{mathpazo}\begin{document} etex;
  
beginfig(7);
  path bobine, bobina, Right, Up, Left, Down;
  numeric scal, outeray, inneray, rotange, rotanga, tens, dist;
  pen grossa;
  scal = 15;
  inneray = 0.5;
  outeray = 2.1;
  rotange = 0;
  rotanga = 12;
  tens = 3.6;
  grossa = pencircle scaled 1pt;
  dist = scal*0.35*(inneray+outeray);
  Right = (inneray*scal*right)--(outeray*scal*right);
  Up = (outeray*scal*up)--(inneray*scal*up);
  Left = (inneray*scal*left)--(outeray*scal*left);
  Down = (outeray*scal*down)--(inneray*scal*down);
  bobine = (Right..tension tens..Up..Left..tension tens..Down..cycle)
  rotated rotange;
  bobina = bobine rotated rotanga;
  draw Right;
  draw Up;
  draw Left;
  draw Down;
  draw fullcircle scaled (2*inneray*scal);
  draw fullcircle scaled (2*outeray*scal);
%  draw bobine dashed evenly withpen grossa;
  draw bobina withpen grossa;
  z1 = dist*(right+up);
  z2 = z1 rotated 90;
  z3 = z1 rotated 180;
  z4 = z1 rotated -90;
  label( btex N etex, z1 );
  label( btex S etex, z2 );
  label( btex N etex, z3 );
  label( btex S etex, z4 );
endfig;

prologues := 1;

beginfig(1);
  numeric u, m, d;
  u = 7.5pt;
  m = 0.09;
  d = 0.25;
  f := (3,5,4);
  Spread := 105;
  TDAtiplen := 0.1;      
  TDAhalftipbase := 0.04; 
  TDAhalfthick := 0.02;
  color v[];
  v1 := (0,0,0);
  v2 := (1,0,0);
  v3 := (1,0,1);
  v4 := (1,1,1);
  v5 := (1,1,0);
  v6 := (0,1,0);
  v7 := (0,1,1);
  v8 := (0,0,1);
  v22 := (1+m,-m,d);
  v23 := (1+m,-m,1-d);
  v25 := (1-d,1+m,-m);
  v26 := (d,1+m,-m);
  v27 := (-m,1-d,1+m);
  v28 := (-m,d,1+m);
  path espiracubica;
  espiracubica = rp(v1)--rp(v2)--rp(v3)--rp(v4)--
    rp(v5)--rp(v6)--rp(v7)--rp(v8)--cycle;
  draw espiracubica withpen pencircle scaled u;
  undraw rp(v4)--rp(v5) withpen pencircle scaled 1.9u;
  draw rp(v3)--rp(v4)--rp(v5)--rp(v6) withpen pencircle scaled u;
  draw rp(v2)--rp(v5) dashed evenly;
  draw rp(v4)--rp(v7) dashed evenly;
  draw rp(v8)--rp(v3) dashed evenly;
  draw rp(v6)--rp(v1) dashed evenly;
  tdarrow(v22,v23);
  tdarrow(v25,v26);
  tdarrow(v27,v28);
endfig;

beginfig(2);
  espiracubica := rp(v1)--rp(v2)--rp(v3)--rp(v4)--rp(v7)--rp(v6)--cycle;
  draw espiracubica withpen pencircle scaled u;
  draw rp(v8)--rp(v1) dashed evenly;
  draw rp(v8)--rp(v3) dashed evenly;
  draw rp(v8)--rp(v7) dashed evenly;
  draw rp(v5)--rp(v2) dashed evenly;
  draw rp(v5)--rp(v4) dashed evenly;
  draw rp(v5)--rp(v6) dashed evenly;
  tdarrow(v22,v23);
endfig;

beginfig(3);
  espiracubica := rp(v1)--rp(v3)--rp(v5)--rp(v7)--cycle;
  v31 := (d-m,0,d+m);
  v33 := (1-d-m,0,1-d+m);
  draw espiracubica withpen pencircle scaled u;
  draw rp(v8)--rp(v1) dashed evenly;
  draw rp(v8)--rp(v3) dashed evenly;
  draw rp(v8)--rp(v7) dashed evenly;
  draw rp(v5)--rp(v2) dashed evenly;
  draw rp(v5)--rp(v4) dashed evenly;
  draw rp(v5)--rp(v6) dashed evenly;
  draw rp(v1)--rp(v2)--rp(v3)--rp(v4)--rp(v7)--rp(v6)--cycle dashed evenly;
  tdarrow(v31,v33);
endfig;

beginfig(4);
  espiracubica := rp(v1)--rp(v8)--rp(v3)--rp(v2)--rp(v5)--rp(v6)--cycle;
  draw espiracubica withpen pencircle scaled u;
  draw rp(v8)--rp(v1) dashed evenly;
  draw rp(v8)--rp(v3) dashed evenly;
  draw rp(v8)--rp(v7) dashed evenly;
  draw rp(v5)--rp(v2) dashed evenly;
  draw rp(v5)--rp(v4) dashed evenly;
  draw rp(v5)--rp(v6) dashed evenly;
  draw rp(v1)--rp(v2)--rp(v3)--rp(v4)--rp(v7)--rp(v6)--cycle dashed evenly;
  tdarrow(v22,v23);
  tdarrow(v26,v25);
endfig;

beginfig(5);
  f := (5,3,4);
  espiracubica := rp(v1)--rp(v2)--rp(v5)--rp(v3)--
    rp(v8)--rp(v7)--rp(v4)--rp(v6)--cycle;
  draw espiracubica withpen pencircle scaled u;
  undraw rp(v3)--rp(v5) withpen pencircle scaled 1.9u;
  draw rp(v2)--rp(v5)--rp(v3)--rp(v8) withpen pencircle scaled u;
  draw rp(v2)--rp(v3)--rp(v4)--rp(v5)--rp(v6)--rp(v7) dashed evenly;
  draw rp(v1)--rp(v8) dashed evenly;
  tdarrow(v27,v28);
endfig;

beginfig(6);
  espiracubica := rp(v1)--rp(v8)--rp(v2)--rp(v3)--
    rp(v5)--rp(v4)--rp(v6)--rp(v7)--cycle;
  draw espiracubica withpen pencircle scaled u;
  undraw rp(v3)--rp(v5) withpen pencircle scaled 1.9u;
  draw rp(v2)--rp(v3)--rp(v5)--rp(v4) withpen pencircle scaled u;
  draw rp(v1)--rp(v2)--rp(v5)--rp(v6)--cycle dashed evenly;
  draw rp(v3)--rp(v4)--rp(v7)--rp(v8)--cycle dashed evenly;
  tdarrow(v22,v23);
endfig;

beginfig(8);
  espiracubica := rp(v2)--rp(v3)--rp(v8)--rp(v7)--
    rp(v5)--cycle;
  draw espiracubica withpen pencircle scaled u;
  draw rp(v3)--rp(v4)--rp(v7)--rp(v6)--rp(v1)--rp(v8) dashed evenly;
  draw rp(v1)--rp(v2) dashed evenly;
  draw rp(v4)--rp(v5)--rp(v6) dashed evenly;
  tdarrow(v23,v22);
  tdarrow(v27,v28);
endfig;

verbatimtex \end{document} etex
  
end.