% geom_freder.mp
% L. Nobre G. 
% 2000

input featpost3Dplus2D;

verbatimtex \documentclass{article} \usepackage{beton} \usepackage{concmath} \usepackage{ccfonts} \begin{document} etex

beginfig(1);

    f := (1.7,1.2,1);

    Spread := 250;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, long, base, thetamax;
    path closingface;
    sizevec := (11,9,7);
    gridvec := (0.1,0.1,0.1);
    thetamax := 50;
    long := 0.08;
    base := 0.1;
    starpos := ( -0.5*X(gridvec)*(X(sizevec)-1) ,
                 -0.5*Y(gridvec)*(Y(sizevec)-1) ,
                 -0.5*Z(gridvec)*(Z(sizevec)-1) );
    bordpos := starpos - 0.5*gridvec;
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
    draw facepath(1);
    counter := 1;
    for i= X(sizevec) upto X(sizevec):
        for j=1 upto Y(sizevec):
            for k=1 upto Z(sizevec):
                phi := 0;
                theta := 0;
                currpos := starpos+( (i-1)*X(gridvec) ,
                                     (j-1)*Y(gridvec) ,
                                     (k-1)*Z(gridvec) );
                phi := thetamax*cosd(90*Z(currpos)/abs(Z(starpos)));
                theta := 90+thetamax*cosd(90*Z(currpos)/abs(Z(starpos)));
                generatedirline( counter, phi, theta, long, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    for i=1 upto X(sizevec)-1:
        for j=Y(sizevec) upto Y(sizevec):
            for k=1 upto Z(sizevec):
                phi := 0;
                theta := 0;
                currpos := starpos+( (i-1)*X(gridvec) ,
                                     (j-1)*Y(gridvec) ,
                                     (k-1)*Z(gridvec) );
                phi := thetamax*cosd(90*Z(currpos)/abs(Z(starpos)));
                theta := 90+thetamax*cosd(90*Z(currpos)/abs(Z(starpos)));
                generatedirline( counter, phi, theta, long, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    NL := decr(counter);
    director_invisible( true, base, false );
    closingface = facepath(2);
    fill closingface withcolor background;
    draw closingface;
    pickup pencircle scaled 1pt;
    currpos := (-X(bordpos)+X(gridvec),Y(bordpos), Z(bordpos));
    starpos := (-X(bordpos)+X(gridvec),Y(bordpos)+4*Y(gridvec), Z(bordpos));
    drawarrow rp(currpos)..rp(starpos);
    label.bot(btex $y$ etex,rp(starpos));
    currpos := (-X(bordpos)+X(gridvec),Y(bordpos), Z(bordpos));
    starpos := (-X(bordpos)+X(gridvec),Y(bordpos), Z(bordpos)+4*Z(gridvec));
    drawarrow rp(currpos)..rp(starpos);
    label.lft(btex $z$ etex,rp(starpos));
    currpos := (X(bordpos),-Y(bordpos)+Y(gridvec), Z(bordpos));
    starpos := (X(bordpos)+4*X(gridvec),-Y(bordpos)+Y(gridvec), Z(bordpos));
    drawarrow rp(currpos)..rp(starpos);
    label.bot(btex $x$ etex,rp(starpos));

endfig;

verbatimtex \end{document} etex

end;

