% verygoodcone.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

beginfig(1);
    numeric radiusin, fac, i, numc, stepang, xcoor, ycoor, radc;
    
    f 			:= (5,0.3,1.5);
    Spread 		:= 70;
    radiusin		= 0.3;
    numc		= 9;
    fac                 = 3;
    stepang             = 360/numc;
    radc                = (fac-1)*radiusin*sind(0.5*stepang);

    for i=0 step stepang until 180:
      xcoor := -radiusin*cosd(i);
      ycoor := radiusin*sind(i);
verygoodcone(true,fac*(xcoor,-ycoor,0),(-xcoor,ycoor,0),radc,(xcoor,-ycoor,0));
verygoodcone(true,fac*(xcoor,ycoor,0),-(xcoor,ycoor,0),radc,(xcoor,ycoor,0));
    endfor;
endfig;

end;


