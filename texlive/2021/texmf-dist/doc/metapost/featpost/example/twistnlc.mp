% explvis2
% L. Nobre G. 
% 1999

input featpost3Dplus2D;

beginfig(1);

    f := (2,1,1);

    Spread := 250;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, thetamin, 
            thetamax, thetastp, long, base, radius, polar;
    sizevec := (4,36,1);
    gridvec := (0.14,0.1,0.3);
    thetamin := 90;
    thetamax := 150;
    long := 0.08;
    base := 0.02;
    thetastp := (thetamax-thetamin)/(X(sizevec)-1);
    starpos := ( -X(gridvec)*X(sizevec) ,
                 -X(gridvec)*X(sizevec) ,
                 0 );
    bordpos := starpos;
    npf1 := 4;
    F1p1 := bordpos;
    F1p2 := (X(bordpos),-Y(bordpos), Z(bordpos));
    F1p3 := (-X(bordpos),-Y(bordpos), Z(bordpos));
    F1p4 := (-X(bordpos),Y(bordpos), Z(bordpos));
    counter := 2;
    for i=1 upto X(sizevec):
        radius := (i-1)*X(gridvec);
        polar := i*floor(Y(sizevec)/X(sizevec));
        for j=1 upto polar:
            for k=1 upto Z(sizevec):
                phi := 360*(j-1)/polar;
                theta := thetamin+(i-1)*thetastp;
                currpos := radius*(cosd(phi-90),sind(phi-90),0)
                           +(k-1)*(0,0,Z(gridvec));
                generatedirface( counter, phi, theta, long, base, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    NF := decr(counter);
    draw_all_test(false);
%    sharpraytrace;
endfig;

end;


