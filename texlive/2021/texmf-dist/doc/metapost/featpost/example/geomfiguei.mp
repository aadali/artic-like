% geomfiguei.mp
% L. Nobre G. 
% 2002

input featpost3Dplus2D;

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts} 
\begin{document} 
etex

beginfig(1);

    f := 3*(1.2,1.7,0.8);

    Spread := 125;

    color origr, nzero, nothr, ndxyl, vmagf, mfxzp, mfxyp, mfxxp;
    numeric alphan, gamman, thetan, phiang, mfsize, axxc, ayyc, azzc;
    numeric phiarmlen, alphafac;
    axxc = 0.99;
    ayyc = 0.87;      % make this 1.28 when in twist-bend geometry
    azzc = 0.82;
    alphan = 170;     % make this 90 to get the twist-bend geometry
    gamman = 25;
    thetan = 10;
    phiang = 30;
    mfsize = 1.25;
    phiarmlen = 0.75;
    alphafac = 1.4;
    origr = (0,0,0);
    nzero = (cosd(gamman),0,sind(gamman));
    ndxyl = phiarmlen*(cosd(phiang),sind(phiang),0);
    nothr = (cosd(phiang)*cosd(gamman+thetan),
             sind(phiang)*cosd(gamman+thetan),sind(gamman+thetan));
    vmagf = mfsize*(cosd(alphan)*sind(gamman),sind(alphan),
                    -cosd(alphan)*cosd(gamman));
    mfxzp = (X(vmagf),0,Z(vmagf));
    mfxyp = (X(vmagf),Y(vmagf),0);
    mfxxp = (X(vmagf),0,0);
    cartaxes( axxc, ayyc, azzc );
    draw rp(origr)--rp(ndxyl) dashed evenly;
    draw rp(origr)--rp(mfxzp)--rp(vmagf) dashed evenly;
    draw rp(origr)--rp(mfxyp)--rp(vmagf) dashed evenly;
    draw rp(origr)--rp(alphafac*mfxyp) dashed evenly;
    draw rp(mfxzp)--rp(mfxxp)--rp(origr) dashed evenly;
    draw rp(mfxyp)--rp(mfxxp) dashed evenly;
    draw rp((0,0,Z(nothr)))--rp(nothr) dashed evenly;
    angline(nzero,(1,0,0),origr,0.65,btex $\gamma$ etex,lft);
    angline(ndxyl,(1,0,0),origr,0.5,btex $\phi$ etex,llft);
    angline(ndxyl,nothr,origr,0.5,btex $\gamma+\theta$ etex,lft);
    angline(vmagf,(0,0,1),origr,0.5,btex $\beta$ etex,top);     % twist-bend
    angline(mfxyp,(1,0,0),origr,0.8,btex $\alpha$ etex,bot);
    pickup pencircle scaled 1.15pt;
    drawarrow rp(origr)..rp(nzero);
    drawarrow rp(origr)..rp(nothr);
    drawarrow rp(origr)..rp(vmagf);
    label.ulft(btex $\vec{n}_0$ etex,rp(nzero));
    label.ulft(btex $\vec{n}$ etex,rp(nothr));
    label.urt(btex $\vec{H}$ etex,rp(vmagf));
endfig;

verbatimtex 
\end{document} 
etex

end;

