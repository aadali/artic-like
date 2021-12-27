% conegetready.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

beginfig(1);
    numeric radiusin, fac, i, xcoor, ycoor, radc, numc, baz;
    color refvec;
    string rv, mrv, vgc, stra, strb, lrv;
    
    f 			:= (5,0.3,1.5);
    Spread 		:= 70;
    radiusin		= 0.3;
    baz                 = 0.08;
    fac                 = 3;
    numc                = 9;
    radc                = (fac-1)*radiusin*sind(180/numc);

    i = 155;
    xcoor := -radiusin*cosd(i);
    ycoor := radiusin*sind(i);
    refvec := (xcoor,-ycoor,baz);
    rv := cstr( refvec );
    mrv := cstr( -refvec );
    lrv := cstr( (2*fac-1)*refvec );
    vgc :="verygoodcone(true,fac*";
    stra := vgc & rv & "," & mrv & "," & decimal(radc) & "," & rv & ")";
    getready( stra, 0.5*(fac+1)*refvec );
    strb := vgc & rv & "," & rv & "," & decimal(radc) & "," & lrv & ")";
    getready( strb, 0.5*(3*fac-1)*refvec );
    doitnow;
    produce_vga_border;
endfig;

end;
