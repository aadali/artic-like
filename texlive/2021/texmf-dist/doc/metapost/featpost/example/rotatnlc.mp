% rotatnlc.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

beginfig(1);

    f := 2*(2,1,0.55);

    Spread := 450;

    color gridvec, sizevec, currpos, starpos, bordpos;
    numeric i, j, k, counter, theta, phi, thetamin; 
    numeric thetamax, thetastp, long, radius, polar;
    sizevec := (7,42,1);       % (number of circles,
                               % number of vectors in outer circle,
                               % number of levels )
    gridvec := (0.08,0.23,0.3);
    thetamin := 90;
    thetamax := 150;
    long := 0.08;
    thetastp := (thetamax-thetamin)/(X(sizevec)-1);
    starpos := ( -X(gridvec)*X(sizevec) ,
                 -X(gridvec)*X(sizevec) ,
                 0 );
    bordpos := starpos;
    counter := 1;
    for i=1 upto X(sizevec):
        radius := (i-1)*X(gridvec);
        polar := i*floor(Y(sizevec)/X(sizevec));
        for j=1 upto polar:
            for k=1 upto Z(sizevec):
                phi := 360*(j-1)/polar;
                theta := thetamin+(i-1)*thetastp;
                currpos := radius*(cosd(phi-90),sind(phi-90),0)
                           +(k-1)*(0,0,Z(gridvec));
                generatedirline( counter, phi, theta, long, currpos );
                counter := incr(counter);
            endfor;
        endfor;
    endfor;
    NL := decr(counter);
    director_invisible( true, Y(gridvec), false );
endfig;

end;


