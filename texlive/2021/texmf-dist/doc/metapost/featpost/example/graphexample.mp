% graphexample.mp
% L. Nobre G.
% 2003

input featpost3Dplus2D;

%input mp-tool;
 
beginfig(1);
    path thezigzag;
    thezigzag = zigzagfrontier( origin, (200,200), 17, 
			4, 20, 16, 12, blue, background );
    draw thezigzag;
endfig;

end;
