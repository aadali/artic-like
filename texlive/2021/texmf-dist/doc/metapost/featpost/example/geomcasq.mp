% geom_casq.mp
% L. Nobre G. 
% 1999

input featpost3Dplus2D;

%input la3Dmacros;
%input hlr3Dmacros;
%input plan2Dmacros;
% 

verbatimtex 
\documentclass{article} 
\usepackage{beton} 
\usepackage{concmath} 
\usepackage{ccfonts}
\begin{document} 
etex

beginfig(1);

    f := (1.7,1.2,1);

    Spread := 150;

    color gridvec, sizevec, currpos, starpos, bordpos, orig, axxc, ayyc, azzc;
    numeric alphamax, long, base, thetamax;
    axxc = (0.8,0,0);
    ayyc = (0,0.7,0);
    azzc = (0,0,0.4);
    sizevec := (11,9,9);
    gridvec := (0.1,0.1,0.1);
    thetamax := 50;
    alphamax := 70;
    long := 0.65;
    base := 0.3;
    starpos := ( -0.5*X(gridvec)*(X(sizevec)-1) ,
                 -0.5*Y(gridvec)*(Y(sizevec)-1) ,
                 -0.5*Z(gridvec)*(Z(sizevec)-1) );
    bordpos := starpos - 0.5*gridvec;
    orig := starpos + 0.5*gridvec + (0.05,0,-0.1);
    npf1 := 4;
    F1p1 := bordpos;
    F1p2 := (X(bordpos),-Y(bordpos), Z(bordpos));
    F1p3 := (-X(bordpos),-Y(bordpos), Z(bordpos));
    F1p4 := (-X(bordpos),Y(bordpos), Z(bordpos));
    npf2 := 4;
    F2p1 := (X(bordpos),Y(bordpos), -Z(bordpos));
    F2p2 := (X(bordpos),-Y(bordpos), -Z(bordpos));
    F2p3 := (-X(bordpos),-Y(bordpos), -Z(bordpos));
    F2p4 := (-X(bordpos),Y(bordpos), -Z(bordpos));
    NF := 2;

    pickup pencircle scaled 1pt;
    drawarrow rp(orig)..rp(orig+axxc);
    drawarrow rp(orig)..rp(orig+ayyc);
    drawarrow rp(orig)..rp(orig+azzc);
    label.bot(btex z etex,rp(orig+axxc));
    label.bot(btex x etex,rp(orig+ayyc));
    label.lft(btex y etex,rp(orig+azzc));
    pickup pencircle scaled 1.5pt;

    starpos := (base,0,0)+orig;
    drawarrow rp(orig)..rp(starpos);
    label.bot(btex $\vec{n}_0$ etex,rp(starpos));

    currpos := starpos;
    starpos := long*(cosd(alphamax),sind(alphamax),0)+orig;
    drawarrow rp(orig)..rp(starpos);
    label.bot(btex $\vec{H}$ etex,rp(starpos));
    angline(starpos,currpos,orig,0.55,btex $\alpha$ etex,bot);

    starpos := base*(cosd(thetamax),sind(thetamax),0)+orig;
    drawarrow rp(orig)..rp(starpos);
    label.bot(btex $\vec{n}$ etex,rp(starpos));
    angline(starpos,currpos,orig,0.2,btex $\theta$ etex,bot);
endfig;

verbatimtex 
\end{document} 
etex

end;

