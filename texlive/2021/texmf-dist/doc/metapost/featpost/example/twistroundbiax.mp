% twistroundbiax.mp
% L. Nobre G. 
% 2002

input featpost3Dplus2D;

 

beginfig(1);

    f := (2,1,1);

    Spread := 250;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, gamma, thetamin, innercircle,
            thetamax, thetastp, long, base, radius, polar, larg;
    sizevec := (4,36,1);
    gridvec := (0.143,0.1,0.1);
    thetamin := 90;
    thetamax := 150;
    long := 0.08;
    larg := 0.02;
    base := 0.2;
    gamma := -30;
    innercircle := 0.5;
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
    NF := 1;
    draw facepath(1);
    counter := 0;
    for i=1 upto X(sizevec):
        radius := (i-innercircle)*X(gridvec);
        polar := i*floor(Y(sizevec)/X(sizevec));
        for j=1 upto polar:
            for k=1 upto Z(sizevec):
                phi := 360*(j-1)/polar;
                theta := thetamin+(i-innercircle)*thetastp;
                currpos := radius*(cosd(phi-90),sind(phi-90),0)
                           +(k-1)*(0,0,Z(gridvec));
                counter := incr(counter);
%                generateonebiax( counter, phi, theta, long, 
%                                          gamma, larg, currpos );
                generatedirline( counter, phi, theta, long, currpos );
            endfor;
        endfor;
    endfor;
    NL := counter;
    director_invisible( true, base, false );
endfig;

end;


