% rigorcubo.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

beginfig(0); 	
	defaultscale := 3;
	defaultfont := "ccr10";
    Spread := 50;
    V1 := (1,1,1);
    V2 := (-1,1,1);
    V3 := (-1,-1,1);
    V4 := (1,-1,1);
    V5 := (1,1,-1);
    V6 := (-1,1,-1);
    V7 := (-1,-1,-1);
    V8 := (1,-1,-1);
    makeface1(1,2,3,4);
    makeface2(5,6,7,8);
    makeface3(1,2,6,5);
    makeface4(2,3,7,6);
    makeface5(3,4,8,7);
    makeface6(4,1,5,8);
    f := (3,5,3);
    pickup pencircle scaled 1pt;
    sharpraytrace;
    label.lft("A", rp(V8) );
    label.top("B", rp(V1) );
    label.top("C", rp(V3) );
    label.lft("D", rp(V4) );
    label.rt("E", rp(V2) );
    label.rt("F", rp(V6) );
endfig;

end;
