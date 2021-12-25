% downloadicon.mp
% L. Nobre G. 
% 2004

input featpost3Dplus2D;

f := 0.4[(12,3,2),(2,5,6)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Max: 0.8
Spread := 30;
TDAtiplen := 0.65;          
TDAhalftipbase := 0.35;     
TDAhalfthick := 0.15;       

beginfig(1);
    color tracecol, lockcol, arrowcol;
    numeric altura, compri, radius, altfrac, comfrac, radang;
    path theside, thefront, thelocka, thelockb, thelock;
    color poi[], startarrow, endarrow;
    numeric numpl, aux, i;
    altura 		= 1;
    compri 		= 2;
    radius		= 0.8;
    altfrac		= 0.2;
    comfrac		= 0.15;
    radang 		= 37;
    numpl 		= 9;
    tracecol 		= blue;
    lockcol 		= black;
    arrowcol            = (0.7,0.2,0.1);
    startarrow          = (0,-radius,2);
    endarrow            = (0,0,radius);
    currentpen := pencircle scaled 2.4mm;
    drawoptions( withcolor tracecol );
    rigorousdisc( 0, true, (-0.5compri,0,0), radius, (compri,0,0) );
    poi1 = (0.5compri,-radius,0);
    poi2 = (0.5compri,radius,0);
    poi3 = (-0.5compri,radius,0);
    poi4 = (0.5compri,-radius,-altura);
    poi5 = (0.5compri,radius,-altura);
    poi6 = (-0.5compri,radius,-altura);
    theside = rp(poi1)--rp(poi2)--rp(poi5)--rp(poi4)--cycle;
    thefront= rp(poi2)--rp(poi3)--rp(poi6)--rp(poi5)--cycle;
    drawoptions( withcolor white );
    fill theside;
    fill thefront;
    thelocka = rp((comfrac*compri,radius,-altfrac*altura))
               --rp((comfrac*compri,radius,0));
    thelockb = rp((-comfrac*compri,radius,-altfrac*altura))
               --rp((-comfrac*compri,radius,0));
    for i=1 upto numpl:
      aux := i*radang/numpl;
      poi7 := (comfrac*compri,radius*cosd(aux),radius*sind(aux));
      poi8 := (-comfrac*compri,radius*cosd(aux),radius*sind(aux));
      thelocka := thelocka...rp(poi7);
      thelockb := thelockb...rp(poi8);
    endfor;
    thelock = thelocka--reverse thelockb--cycle;
    drawoptions( withcolor lockcol );
    filldraw thelock;
    drawoptions( withcolor tracecol );
    draw theside;
    draw thefront;
    drawoptions( withcolor arrowcol );
    tdarrow( startarrow, endarrow );
endfig;

end;
