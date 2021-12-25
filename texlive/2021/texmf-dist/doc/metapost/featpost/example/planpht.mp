% planpht.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

% Planification wise 3D figures

verbatimtex
\documentclass{article}
\usepackage{beton}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

beginfig(1); % phtpbq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definition of the 3D-figure
    V1 := (1,0,0);
    V2 := (0,1,0);
    V5 := (0,0,1);
    V3 := (-1,0,0);
    V4 := (0,-1,0);
    makeface1(1,2,3,4);
    makeface2(1,2,5);
    makeface3(2,3,5);
    makeface4(3,4,5);
    makeface5(4,1,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := (3,5,4);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
endfig;    

beginfig(2); % phtdpbp
    numeric altura, raio, xp, xs, yp, ys;
    altura = 0.5257311;
    raio = 0.8506508;
    yp = 0.80901699;
    ys = 0.5;
    xp = 0.2628655;
    xs = -0.6881909;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definition of the 3D-figure
    V1 := (raio,0,0);
    V2 := (xp,yp,0);
    V3 := (xs,ys,0);
    V4 := (xs,-ys,0);
    V5 := (xp,-yp,0);
    V6 := (0,0,altura);
    V7 := (0,0,-altura);
    makeface1(1,2,6);
    makeface2(5,1,6);
    makeface3(2,3,6);
    makeface4(3,4,6);
    makeface5(4,5,6);
    makeface6(1,2,7);
    makeface7(5,1,7);
    makeface8(2,3,7);
    makeface9(3,4,7);
    makeface10(4,5,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := (3,5,1.8);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
endfig;    

beginfig(3); % phtcubo eps
	defaultscale := 1;
	defaultfont := "ccr10";
    V1 := (1,1,1);
    V2 := (-1,1,1);
    V3 := (-1,-1,1);
    V4 := (1,-1,1);
    V5 := (1,1,-1);
    V6 := (-1,1,-1);
    V7 := (-1,-1,-1);
    V8 := (1,-1,-1);
    makeface1(1,2,3,4);
    makeface2(5,6,7,8);
    makeface3(1,2,6,5);
    makeface4(2,3,7,6);
    makeface5(3,4,8,7);
    makeface6(4,1,5,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := (3,5,3);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
    label.lft("A", rp(V8) );
    label.top("B", rp(V1) );
    label.top("C", rp(V3) );
    label.lft("D", rp(V4) );
    label.rt("E", rp(V2) );
    label.rt("F", rp(V6) );
endfig;

beginfig(4); % phtptq
    numeric high, low;
    high = sqrt(2);
    low = sqrt(2)-1.6818;
    V1 := (1,1,high);
    V2 := (-1,1,high);
    V3 := (-1,-1,high);
    V4 := (1,-1,high);
    V5 := (high,0,low);
    V6 := (0,high,low);
    V7 := (-high,0,low);
    V8 := (0,-high,low);
    V9 := (0,0,0);
    makeface1(5,6,7,8);
    makeface2(6,1,2);
    makeface3(7,2,3);
    makeface4(8,3,4);
    makeface5(5,4,1);
    makeface6(5,6,1);
    makeface7(6,7,2);
    makeface8(7,8,3);
    makeface9(8,5,4);
    makeface10(9,1,2);
    makeface11(9,2,3);
    makeface12(9,3,4);
    makeface13(9,4,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := (2,6,5);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
endfig;

beginfig(5); % phtnonrec
    numeric tro;
    tro = 2;
    V1 := (tro,0,0);
    V2 := (0,tro,0);
    V3 := (-tro,0,0);
    V4 := (0,-tro,0);
    V5 := (1,1,1);
    V6 := (-1,1,1);
    V7 := (-1,-1,1);
    V8 := (1,-1,1);
    V9 := (1,1,-1);
    V10 := (-1,1,-1);
    V11 := (-1,-1,-1);
    V12 := (1,-1,-1);
    V13 := (0,0,tro);
    V14 := (0,0,-tro);
    makeface1(1,5,13,8);
    makeface2(2,6,13,5);
    makeface3(3,7,13,6);
    makeface4(4,8,13,7);
    makeface5(1,5,2,9);
    makeface6(2,6,3,10);
    makeface7(3,7,4,11);
    makeface8(4,8,1,12);
    makeface9(1,9,14,12);
    makeface10(2,10,14,9);
    makeface11(3,11,14,10);
    makeface12(4,12,14,11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f := 2*(5,3,1.5);
    Spread := 50;
    pickup pencircle scaled 1pt;
    sharpraytrace;
endfig;

beginfig(6); % pauphtflex eps
	defaultscale := 1;
	defaultfont := "ccr10";
    color A, B, C, D, E, F, G;
    numeric tro, Adist;
    tro = sqrt(3)/2;
    Adist = 0.88854;
    V1 := (-1.8,0,0);
    V2 := (-1.8,0.6,1.2*tro);
    V3 := (1,0.6,1.2*tro);
    V4 := (1,0,0);
    V5 := (1,2.5,0);
    V6 := (-1.8,2.5,0);
    V7 := (0,0,0);
    V8 := (0,0.707,0);
    V9 := (0,0.707/4,0.707*tro/2);
    V10:= (0.21928,Adist/2,tro*Adist);
    V11:= (-0.37392,0.39638,0);
    V12:= (-1.2125,Adist,0);
    makeface1(1,2,3,4);
    makeface2(4,5,6,1);
    makeface3(7,8,9);
    makeface4(8,9,10);
    makeface5(8,10,11);
    makeface6(10,11,12);
    A = (0,1.3,0);
    B = (0,1.3+sqrt(0.75),0);
    C = (-0.5,1.3+sqrt(0.75),0);
    D = (E+B)/2;
    E = (sqrt(0.75),1.3,0);
    F = ((1+sqrt(3))/4,1.3+(1+sqrt(3))/4,0);
    G = (-1.5,1.3+sqrt(0.75),0);
    f := 3*(2,5,2);
    Spread := 100;
    draw rp(A)--rp(B)--rp(E)--cycle;
    draw rp(F)--rp(B)--rp(E)--cycle;
    draw rp(A)--rp(B)--rp(C)--cycle;
    draw rp(A)--rp(C)--rp(G)--cycle;
    draw rp(A)--rp(F) dashed evenly;
    pickup pencircle scaled 1pt;
    sharpraytrace;
    draw rp(V7)--rp(V9);
    draw rp(V10)--rp(V9);
    draw rp(V12)--rp((-1.2125,0,0)) dashed evenly;
    label.urt("Adist", rp((V12+(-1.2125,0,0))/2) );
    label.lft("A", rp(V10) );
    label.bot("B", rp(V8) );
    label.bot("B", rp(B) );
    label.lrt("C", rp(V11) );
    label.lft("D", rp(V9) );
    label.bot("D", rp(D) );
    label.llft("F", rp(V7) );
    label.bot("F", rp(F) );
    label.rt("G", rp(V12) );
    angline(V6, V2, V1, 0.6, btex $\alpha$ etex, rt);
    angline(V10, 2*V9, V9, 0.7, btex $\theta$ etex, bot);
    draw rp(V9+N(V10-V9)*0.65)--rp(V9+N(V10-V9)*0.8);
    draw rp(V9+N(V9)*0.5)--rp(V9+N(V9)*0.8);
    label.bot(btex $p$ etex, rp((F+B)/2) );
    label.ulft(btex $q$ etex, rp((D+B)/2) );
    label.lft(btex $r$ etex, rp((D+A)/2) );
    label.llft(btex $s$ etex, rp((A+B)/2) );
    label.llft(btex $t$ etex, rp((A+C)/2) );
    label.bot(btex $t$ etex, rp((G+C)/2) );
    label.bot(btex $u$ etex, rp((C+B)/2) );
    label.top(btex $v$ etex, rp((A+G)/2) );
endfig;

beginfig(7); % pauphttoro eps
	defaultscale := 1;
	defaultfont := "ccr10";
    numeric sinthe, costhe, ux, uz, lx, lz, co;
    costhe = (4*sqrt(5)-sqrt(3))/11;
    sinthe = sqrt(1-costhe**2);
    ux = costhe*sqrt(2)/2;
    uz = sinthe/2;
    lx = sqrt(4/3)*sinthe;
    lz = -sqrt(2/3)*costhe;
    co = sqrt(2)*(2*sqrt(5)+5*sqrt(3))/22;
    V1 := (ux,0,-uz);
    V2 := (lx,0,-lz);
    V3 := (0,ux,uz);
    V4 := (0,lx,lz);
    V5 := (co,co,0);
    V6 := (-ux,0,-uz);
    V7 := (-lx,0,-lz);
    V8 := (-co,co,0);
    makeface1(1,2,3);
    makeface2(1,2,5);
    makeface3(2,3,5);
    makeface4(3,4,5);
    makeface5(3,4,1);
    makeface6(1,4,5);
    makeface7(3,6,7);
    makeface8(3,7,8);
    makeface9(6,7,8);
    makeface10(3,8,4);
    makeface11(4,6,8);
    makeface12(3,4,6);
    Spread := 100;
%    f := (1,-0.4,3);
    f := 2*(2,1,2.8);
    draw_invisible( false, false, TableC6, black );
%    cartaxes(0.2,0.2,0.2);
    draw rp(V2)--rp(V1)--rp(V4) dashed evenly;
    draw rp(V3)--rp(V1)--rp(V5) dashed evenly;
    draw rp(V6)--rp(V4) dashed evenly;
    draw rp(V6)--rp(V8) dashed evenly;
%    pickup pencircle scaled 1pt;
%    sharpraytrace;
    label.rt("A", rp(V5) );
    label.rt("C", rp(V3+(0,0.05,0)) );
    label.rt("B", rp(V4) );
    label.rt("D", rp(V8) );
    label.lft("E", rp(V2) );
    label.llft("F", rp(V1) );
    label.lft("G", rp(V7) );
    label.lft("H", rp(V6) );
endfig;

verbatimtex \end{document} etex

end;
 
