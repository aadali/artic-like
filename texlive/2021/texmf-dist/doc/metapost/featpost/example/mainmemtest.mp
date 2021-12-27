% mainmemtest.mp
% L. Nobre G. 
% 2007

input featpost3Dplus2D;

beginfig(1);

    f := (2,1,1);

    Spread := 500;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, thetamin, 
            thetamax, thetastp, long, base;
% It should work with sizevec := (8,6,1);
% with (9,6,1) and sharpraytrace it says
% ! MetaPost capacity exceeded, sorry [main memory size=2199001].
%    sizevec := (7,6,1);
% This has been degrading all along the growth of FeatPost.
% Now the message is
% ! MetaPost capacity exceeded, sorry [independent variables=33554431].
%    sizevec := (6,6,1);
    sizevec := (4,4,1);
    gridvec := (0.05,0.1,0.1);
    thetamin := 55;
    thetamax := 125;
    long := 0.08;
    base := 0.01;
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
    counter := 2;
    for i=1 upto X(sizevec):
        for j=1 upto Y(sizevec):
            for k=1 upto Z(sizevec):
                phi := 90;
                theta := thetamin+(i-1)*thetastp;
                currpos := starpos+( (i-1)*X(gridvec) ,
                                     (j-1)*Y(gridvec) ,
                                     (k-1)*Z(gridvec) );
                generatedirface( counter, phi, theta, long, base, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    NF := decr(counter);
    sharpraytrace;
endfig;

end;

