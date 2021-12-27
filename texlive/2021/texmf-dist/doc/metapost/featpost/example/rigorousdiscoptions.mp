% rigorousdiscoptions.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

f := 0.8*(6,3,4);

beginfig(1);
    color axe, cstep, colen;
    numeric nstep, sfrac, ouray, inray;
    
    nstep = 0.65;
    sfrac = 0.1;
    ouray = 0.7;
    inray = 0.5;

    axe = red;
    cstep = ( nstep, 0, 0 );
    colen = ( nstep*sfrac, 0, 0 );

    pickup pencircle scaled 2pt;
    rigorousdisc( 0, true, -cstep, ouray, colen );
    rigorousdisc( 0, false, black, ouray, colen );
    rigorousdisc( inray, true, cstep, ouray, colen );
    rigorousdisc( inray, false, 2*cstep, ouray, colen );
endfig;

end;


