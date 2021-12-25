% polyhedr.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{newcent}
\pagestyle{empty}
\begin{document}
etex

% This exemplifies hidden line removal. 

f := (2,5,3);

beginfig(1);
    Print_Step := 1;
    faceraytrace( 0.4, red );
    lineraytrace( 0.4, black );
%    pickup pencircle scaled 2pt;
%    draw_all_test( false );
%    pickup pencircle scaled 1pt;
%    sharpraytrace;
endfig;    

f := 2*(2,5,3.5);
Print_step := 10;
Spread := 156;
LightSource := 2.5*(4,-3,17);

beginfig(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definition of the 3D-figure
    V1 := (1.4,0,0);
    V2 := (1,1,0);
    V3 := (0,1.4,0);
    V4 := (-1,1,0);
    V5 := (-1.4,0,0);
    V6 := (-1,-1,0);
    V7 := (0, -1.4,0);
    V8 := (1,-1,0); 
    V9 := (1.4,0,1);
    V10:= (1,1,1);
    V11:= (0,1.4,1);
    V12:= (-1,1,1);
    V13:= (-1.4,0,1);
    V14:= (-1,-1,1);
    V15:= (0,-1.4,1);
    V16:= (1,-1,1); 

    V17:= (1.4,0,-1);
    V18:= (1,1,-1);
    V19:= (0,1.4,-1);
    V20:= (-1,1,-1);
    V21:= (-1.4,0,-1);
    V22:= (-1,-1,-1);
    V23:= (0,-1.4,-1);
    V24:= (1,-1,-1);
    
    V25:= (0.7,0,0);
    V26:= (0,0.7,0);
    V27:= (0.7,-0.7,0.7); 
    V28:= (0.7,0,0.7);
    V29:= (0,0.7,0.7);
    V30:= (-0.7,0.7,0.7);
    V31:= (0.35,0.35,1.05);
    V32:= (0,-0.7,1.4);
    V34:= (0,0,1.4);
    V33:= (-0.7,0,1.4); 
    NL := 1;
    npl1 := 9; 
    L1p1 := V17;
    L1p2 := V18;
    L1p3 := V19;
    L1p4 := V20;
    L1p5 := V21;
    L1p6 := V22;
    L1p7 := V23;
    L1p8 := V24;
    L1p9 := V17;
    NF := 8;
    npf1 := 4;
    F1p1 := V1;
    F1p2 := V2;
    F1p3 := V10;
    F1p4 := V9;
    npf2 := 4;
    F2p1 := V2;
    F2p2 := V3;
    F2p3 := V11;
    F2p4 := V10;
    npf3 := 4;
    F3p1 := V3;
    F3p2 := V4;
    F3p3 := V12;
    F3p4 := V11;
    npf4 := 4;
    F4p1 := V4;
    F4p2 := V5;
    F4p3 := V13;
    F4p4 := V12;
    npf5 := 4;
    F5p1 := V5;
    F5p2 := V6;
    F5p3 := V14;
    F5p4 := V13;
    npf6 := 4;
    F6p1 := V6;
    F6p2 := V7;
    F6p3 := V15;
    F6p4 := V14;
    npf7 := 4;
    F7p1 := V7;
    F7p2 := V8;
    F7p3 := V16;
    F7p4 := V15;
    npf8 := 4;
    F8p1 := V8;
    F8p2 := V1;
    F8p3 := V9;
    F8p4 := V16;
    makeface9(25,26,29,31,28);
    makeface10(25,27,28);
    makeface11(26,29,30);
    makeface12(27,28,31,34,32);
    makeface13(31,29,30,33,34);
    makeface14(32,33,34);
    makeface15(25,26,30,33,32,27);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%    draw_all_test(true);
    ShadowOn := true;
    fill_faces();
    draw_invisible( true, true, TableC3, TableC0 );
%    lineraytrace(0.3, black);
%    faceraytrace(0.4, green);			% This my take very long.
%    sharpraytrace;

endfig;    

f := 2*(-2,-0.5,2.5);
Spread := 56;
LightSource := 1.5*(-4,5,2.7);

beginfig(3);
  numeric i, aux;
  pair anglepar;
  V1 := ( 0.7, 0,   0   );
  V2 := ( 0,   0.7, 0   );
  V3 := ( 0.7,-0.7, 0.7 ); 
  V4 := ( 0.7, 0,   0.7 );
  V5 := ( 0,   0.7, 0.7 );
  V6 := (-0.7, 0.7, 0.7 );
  V7 := ( 0.35,0.35,1.05);
  V8 := ( 0,  -0.7, 1.4 );
  V9 := ( 0,   0,   1.4 );
  V10:= (-0.7, 0,   1.4 );
  anglepar = getanglepair( (1,1,1) );
  for i=11 upto 20:
    V[i]:=eulerrotation( -45, 0, 0, V[i-10] );
    V[i]:=eulerrotation( 0, 90-ypart anglepar, 0, V[i] );
  endfor;
  aux = Z( V11 );
  for i=25 upto 34:
    V[i]:=( X(V[i-14]),Y(V[i-14]),1.35*(Z(V[i-14])-aux));
  endfor;
  makeface1(25,26,29,31,28);
  makeface2(25,27,28);
  makeface3(26,29,30);
  makeface4(32,33,34);
  makeface5(25,26,30,34,32,27);
  makeface6(27,28,31,33,32);
  makeface7(31,29,30,34,33);
  for i=1 upto 7:
    FC[i] := 4;
    FCD[i]:= true;
  endfor;
  ShadowOn := true;
  setthearena( 10, 5 );
  fill_faces();
  draw_invisible( false, false, TableC4, TableC0 );
%   for i=25 upto 34:
%     label.top( decimal(i), rp(V[i]) );
%   endfor;
endfig;    

verbatimtex \end{document} etex

end;
 
