% cubicfacecentered.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

% verbatimtex 
% \documentclass{article} 
% \usepackage{beton} 
% \usepackage{concmath} 
% \usepackage{ccfonts} 
% \begin{document} 
% etex

beginfig(1); 
    numeric tro;
    pen fina, media, grossa;
    color mygray;
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
    V15 := (tro,0,tro);
    V16 := (0,tro,tro);
    V17 := (tro,tro,0);
    V18 := (tro,0,-tro);
    V19 := (tro,-tro,0);
    V20 := (0,-tro,tro);
    V21 := (-tro,0,tro);
    V22 := (-tro,tro,0);
    V23 := (0,tro,-tro);
    V24 := (-tro,0,-tro);
    V25 := (0,-tro,-tro);
    V26 := (-tro,-tro,0);
    V27 := (-2*tro,-tro,-tro);
    V28 := (-2*tro,-tro,tro);
    V29 := (-2*tro,tro,tro);
    V30 := (-2*tro,tro,-tro);
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

    makeface13(27,28,29,30);
    makeface14(20,21,16,29,28);
    makeface15(23,24,25,27,30);
    makeface16(16,22,23,30,29);
    makeface17(25,27,28,20,26);
    makeface18(21,22,24,26);
    makeface19(20,21,26);
    makeface20(16,21,22);
    makeface21(22,23,24);
    makeface22(24,25,26);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    f := 2*(5,3,1.5);
%    f := 2*(5,3,4);
    f := 2.34*(15,2.4,5.05);
    Spread := 50;
    fina  = pencircle scaled 1pt;
    media = pencircle scaled 1.2pt;
    grossa= pencircle scaled 3pt;
    mygray = 0.3background;
%    pickup fina;
%    draw rp(V25)--rp(V27)--rp(V28)--rp(V20)--
%         rp(V25)--rp(V23)--rp(V30)--rp(V27);
    pickup media;
    draw rp(V24)--rp(V25)--rp(V19)--rp(V20)--
    rp(V26)--rp(V25)--rp(V18)--rp(V23)--
    rp(V24)--rp(V26)--rp(V21)--rp(V22)--cycle withcolor mygray;
    drawoptions( withpen grossa );
    draw_invisible( false, false, TableC3, black );
    fillfacewithlight( 14 );
    fillfacewithlight( 16 );
    drawoptions();
    pickup media;
    draw rp(V15)--rp(V16)--rp(V22)--rp(V23)--
         rp(V17)--rp(V16)--rp(V21)--rp(V20)--
         rp(V15)--rp(V17)--rp(V18)--rp(V19)--cycle withcolor mygray;
%    pickup fina;
%    draw rp(V28)--rp(V29)--rp(V30);
%    draw rp(V20)--rp(V16)--rp(V23);
%    draw rp(V16)--rp(V29);

endfig;

% verbatimtex \end{document} etex

end;





