% twistflat.mp
% L. Nobre G. 
% Y2K

input featpost3Dplus2D;

 

beginfig(1);

    f := 0.7*(1.5,0.3,1.2);

    Spread := 700;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, thetamin, 
            thetamax, thetastp, long, base;
    sizevec := (11,6,1);
    gridvec := (0.05,0.1,0.1);
    thetamin := 55;
    thetamax := 125;
    long := 0.08;
    base := 0.1;
    thetastp := (thetamax-thetamin)/(X(sizevec)-1);
    starpos := ( -0.5*X(gridvec)*(X(sizevec)-1) ,
                 -0.5*Y(gridvec)*(Y(sizevec)-1) ,
                 -0.5*Z(gridvec)*(Z(sizevec)-1) );
    bordpos := starpos - 0.5*gridvec;
    npf1 := 4;
    F1p1 := bordpos;
    F1p2 := (X(bordpos),-Y(bordpos), Z(bordpos));
    F1p3 := (-X(bordpos),-Y(bordpos), Z(bordpos));
    F1p4 := (-X(bordpos),Y(bordpos), Z(bordpos));
    NF := 1;
    draw facepath(1);
    counter := 1;
    for i=1 upto X(sizevec):
        for j=1 upto Y(sizevec):
            for k=1 upto Z(sizevec):
                phi := 90;
                theta := thetamin+(i-1)*thetastp;
                currpos := starpos+( (i-1)*X(gridvec) ,
                                     (j-1)*Y(gridvec) ,
                                     (k-1)*Z(gridvec) );
                generatedirline( counter, phi, theta, long, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    NL := decr(counter);
    director_invisible( true, base, false );
endfig;

end;

