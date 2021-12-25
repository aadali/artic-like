% geommine.mp
% L. Nobre G. 
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
%\usepackage{beton} 
%\usepackage{concmath} 
%\usepackage{ccfonts} 
%\usepackage{bbold}
\begin{document} 
etex

prologues := 3;

beginfig(1);

    beta := 15;
    alpha := 40;
    psi := 26;
    theta := 32;
    phi := -33;
    h := 1.35;
    v := 0.4;
    xlen := 1.2;
    ylen := 0.65;
    phirad := 0.67;

    f := (3,5,3);
    Spread := 100;

    color c[], n[], H, h[], Vel, xloc, yloc, vx, vy;
    c0 := (0,0,0);
    cartaxes(xlen,ylen,1.0);
    n0 := (cosd(beta),0,sind(beta));
    n3 := (0,sind(phi),0);
    n2 := n0 + n3;
    n5 := n3 + cosd(phi)*n0;
    n4 := cosd(theta)*n5;
    Fx:=sind(alpha)*cosd(beta)-cosd(alpha)*cosd(psi)*sind(beta);
    Fy:=cosd(alpha)*sind(psi);
    Fz:=sind(alpha)*sind(beta)+cosd(alpha)*cosd(psi)*cosd(beta);
    nx:=cosd(theta)*cosd(phi)*cosd(beta)-sind(theta)*sind(beta);
    ny:=cosd(theta)*sind(phi);
    nz:=cosd(theta)*cosd(phi)*sind(beta)+sind(theta)*cosd(beta);
    n1 := (nx,ny,nz);
    xloc := xlen*n1;
    yloc := ylen*(-greenpart n1,redpart n1,0); 
    H := h*(Fx,Fy,Fz);
    Vel := v*(Fx,Fy,Fz);
    vx := v*sind(phi)*n1;
    vy := v*cosd(phi)*N(yloc);
    h4 := h*sind(phi)*n1;
    h5 := h*cosd(phi)*N(yloc);
    h2 := (-h*sind(beta)*cosd(alpha)*cosd(psi),
           0,
           h*cosd(beta)*cosd(alpha)*cosd(psi));
    h3 := (0,h*cosd(alpha)*sind(psi),0);
    h1 := h2 + h3;
    draw rp(n0)--rp(n2)--rp(n3)--rp(c0);
    draw rp(c0)--rp(n4)--rp(n1);
    draw rp(n5)--rp(n4);
    draw rp(c0)--rp(h2)--rp(h1)--rp(h3)--cycle;
    draw rp(c0)--rp(h1)--rp(H);
%    draw rp(h4)--rp(H)--rp(h5);
    label.llft(btex $\vec{n}_0$ etex,rp(n0));
    label.urt(btex $\vec{n}$ etex,rp(n1));
    label.llft(btex $\vec{H}$ etex,rp(H));
%    label.urt(btex $\vec{v}$ etex,rp(Vel));
%    label.bot(btex $\mathbb{x}$ etex,rp(xloc));
%    label.bot(btex $\mathbb{y}$ etex,rp(yloc));
%    label.lrt(btex $\vec{v}_\mathbb{x}$ etex,rp(vx));
%    label.urt(btex $\vec{v}_\mathbb{y}$ etex,rp(vy));
%    label.lft(btex $\vec{H}_\mathbb{x}$ etex,rp(h4));
%    label.urt(btex $\vec{H}_\mathbb{y}$ etex,rp(h5));
%    pickup pencircle scaled 1pt;
    angline(n0,n5,c0,phirad,btex $\phi$ etex,lft);
%    angline((0,1,0),yloc,c0,phirad,btex $\phi$ etex,rt);
    angline(h2,h1,c0,0.5,btex $\psi$ etex,urt);
    angline(H,h1,c0,0.9,btex $\alpha$ etex,top);
    angline(n5,n1,c0,0.6,btex $\theta$ etex,ulft);
    angline((0,0,1),h2,c0,0.5,btex $\beta$ etex,top);
    angline((1,0,0),n0,c0,0.6,btex $\beta$ etex,rt);
%    drawarrow rp(c0)..rp(xloc);
%    drawarrow rp(c0)..rp(yloc);
    pickup pencircle scaled 1.3pt;
    drawarrow rp(c0)..rp(n0);
    drawarrow rp(c0)..rp(n1);
    drawarrow rp(c0)..rp(H);
%    drawarrow rp(c0)..rp(Vel);
%    drawarrow rp(c0)..rp(vx);
%    drawarrow rp(c0)..rp(vy);
%    drawarrow rp(c0)..rp(h4);
%    drawarrow rp(c0)..rp(h5);
endfig;

verbatimtex \end{document} etex

end;





