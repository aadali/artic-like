% geombasic.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

verbatimtex
\documentclass{article}
\usepackage{beton}
\usepackage{concmath}
\usepackage{ccfonts}
\begin{document}
etex

% Kind of graph that physicists do a lot. Play with the angles
% and h a little to see what happens. Change the point of view f.
% Change the position of symbols using the last argument of angline. 

beginfig(1);

    f := (3,5,2);

    Spread := 155; 

    ShiftV := origin;  

    gamma := 30;
    theta := 15;
    phi := 40;
    h := 1.2;

    p := 0.75;
    ppp := 0.8;

    mythick := 0.9pt;

    color c[], n[], H;
    c0 := (0,0,0);
    n0 := (cosd(gamma),0,sind(gamma));
    n1 := (cosd(gamma+theta),0,sind(gamma+theta));
    n2 := (cosd(gamma)*cosd(phi),cosd(gamma)*sind(phi),sind(gamma));
    n3 := (cosd(gamma+theta)*cosd(phi),
           cosd(gamma+theta)*sind(phi),sind(gamma+theta));
    n4 := (cosd(phi),sind(phi),0);
    H := h*(-sind(gamma),0,cosd(gamma));
    c1 := (0,0,Z(p*n1));
    c2 := (0,0,Z(p*n2));
    c3 := (1,0,0);
    c4 := (0,0,1);
    c5 := (X(H),0,0);
    pickup pencircle scaled 1.25pt;
    drawarrow rp(c0)..rp(n0);
    drawarrow rp(c0)..rp(n3);
    drawarrow rp(c0)..rp(H);
    label.urt(btex $\vec{n}_0$ etex,rp(n0));
    label.urt(btex $\vec{n}$ etex,rp(n3));
    label.lft(btex $\vec{H}$ etex,rp(H));
    pickup pencircle scaled mythick;
    cartaxes(1,1,1);
%    draw rp(c0)--rp(p*n1) dashed evenly scaled 2;
%    draw rp(c0)--rp(p*n2) dashed evenly scaled 2;
    draw rp(c0)--rp(c5)--rp(H) dashed evenly scaled 1.5;
%    draw rp(c1)..rp(p*n1) dashed evenly scaled 2;
    draw rp(c1)..rp(p*n3) dashed evenly scaled 1.5;
%    draw rp(c2)..rp(p*n0) dashed evenly scaled 2;
%    draw rp(c2)..rp(p*n2) dashed evenly scaled 2;
    draw rp(c0)..rp(n4) dashed evenly scaled 1.5;
    angline(c3,n0,c0,ppp,btex $\gamma$ etex,lft);
%    angline(n0,n1,c0,p,btex $\theta$ etex,rt);
%    angline(p*n0,p*n2,c2,X(p*n0),btex $\phi$ etex,lrt);
%    angline(p*n1,p*n3,c1,X(p*n1),btex $\phi$ etex,ulft);
%    angline(n2,n3,c0,p,btex $\theta $ etex,urt);
    angline(c4,H,c0,ppp,btex $\gamma$ etex,top);
    angline(n4,n3,c0,p,btex $\gamma+\theta$ etex,lft);
    angline(n4,c3,c0,p,btex $\phi$ etex,llft);
endfig;

verbatimtex
\end{document}
etex

end.
