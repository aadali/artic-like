% srinathpolineni.mp
% L. Nobre G.
% 2013

input featpost3Dplus2D;

numeric thickn;
thickn = 7;

beginfig(8);
  f := 3*(5,3,4);
  Spread := 60;
  
  numeric hstep, hmax, hsdev, basestep, basemax, basesdev;
  numeric i, j, k, angsdev, cylength, cyradius, basen, hn;
  numeric vcx, vcy, vcz, hcurr, xcurr, ycurr, hbase, aone;
  numeric atwo, zcurr, counter, aonebase, atwobase;
  color lenvec, currpos;
  
  cylength   := 0.6;
  cyradius   := 0.1;
  basen	     := 12;
  hn	     := 4;
  basestep   := cyradius*3.1;
  hstep	     := cylength*1.5;
  basesdev   := cyradius*0.1;
  hsdev	     := hstep*0.2;
  hbase	     := -0.8;
  angsdev    := 7;
  aonebase   := 170;
  atwobase   := 12;
  
  basemax    := basen*basestep;
  hmax	     := hn*hstep;
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      hcurr := hbase;
      forever:
	exitunless hcurr < hmax+hbase;
	hcurr   := hcurr + hstep;
	zcurr   :=         hcurr + hsdev*normaldeviate;
	xcurr   := i*basestep + basesdev*normaldeviate;
	ycurr   := j*basestep + basesdev*normaldeviate;
	aone    :=               uniformdeviate( 360 );
	atwo    :=               angsdev*normaldeviate;
	vcz     := cosd( atwo );
	vcy     := sind( atwo )*sind( aone );
	vcx     := sind( atwo )*cosd( aone );
	currpos := ( xcurr, ycurr, zcurr );
	lenvec  := cylength*( vcx, vcy, vcz ); 
	counter := incr( counter );
	generatedirline( counter, aone, 90-atwo, cylength, currpos );
      endfor;
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(5);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      for hcurr=hbase step hstep until hmax:
	zcurr   := hcurr;
	xcurr   := i*basestep;
	ycurr   := j*basestep;
	aone    := aonebase;
	atwo    := atwobase;
	vcz     := cosd( atwo );
	vcy     := sind( atwo )*sind( aone );
	vcx     := sind( atwo )*cosd( aone );
	currpos := ( xcurr, ycurr, zcurr );
	lenvec  := cylength*( vcx, vcy, vcz ); 
	counter := incr( counter );
	generatedirline( counter, aone, 90-atwo, cylength, currpos );
      endfor;
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(6);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      for hcurr=hbase step hstep until hmax:
	zcurr   := hcurr + basesdev*normaldeviate;
	xcurr   := i*basestep + basesdev*normaldeviate;
	ycurr   := j*basestep + basesdev*normaldeviate;
	aone    := aonebase + angsdev*(normaldeviate-0.5);
	atwo    := atwobase + angsdev*(normaldeviate-0.5);
	vcz     := cosd( atwo );
	vcy     := sind( atwo )*sind( aone );
	vcx     := sind( atwo )*cosd( aone );
	currpos := ( xcurr, ycurr, zcurr );
	lenvec  := cylength*( vcx, vcy, vcz ); 
	counter := incr( counter );
	generatedirline( counter, aone, 90-atwo, cylength, currpos );
      endfor;
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(7);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      for hcurr=hbase step hstep until hmax:
	zcurr   := hcurr + basesdev*normaldeviate;
	xcurr   := i*basestep + basesdev*normaldeviate;
	ycurr   := j*basestep + basesdev*normaldeviate;
	aone    := aonebase + angsdev*(normaldeviate-0.5);
	atwo    :=            angsdev*(normaldeviate-0.5);
	vcz     := cosd( atwo );
	vcy     := sind( atwo )*sind( aone );
	vcx     := sind( atwo )*cosd( aone );
	currpos := ( xcurr, ycurr, zcurr );
	lenvec  := cylength*( vcx, vcy, vcz ); 
	counter := incr( counter );
	generatedirline( counter, aone, 90-atwo, cylength, currpos );
      endfor;
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(9);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      for hcurr=hbase step hstep until hmax:
	zcurr   := hbase + uniformdeviate(hmax-hbase);
	xcurr   := uniformdeviate(basen*basestep);
	ycurr   := uniformdeviate(basen*basestep);
	aone    := uniformdeviate(180);
	atwo    := uniformdeviate(180)-90;
	vcz     := cosd( atwo );
	vcy     := sind( atwo )*sind( aone );
	vcx     := sind( atwo )*cosd( aone );
	currpos := ( xcurr, ycurr, zcurr );
	lenvec  := cylength*( vcx, vcy, vcz ); 
	counter := incr( counter );
	generatedirline( counter, aone, 90-atwo, cylength, currpos );
      endfor;
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(1);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      zcurr   := basesdev*normaldeviate;
      xcurr   := uniformdeviate(basen*basestep);
      ycurr   := uniformdeviate(basen*basestep);
      aone    := uniformdeviate(180);
      atwo    := 90;
      vcz     := cosd( atwo );
      vcy     := sind( atwo )*sind( aone );
      vcx     := sind( atwo )*cosd( aone );
      currpos := ( xcurr, ycurr, zcurr );
      lenvec  := cylength*( vcx, vcy, vcz ); 
      counter := incr( counter );
      generatedirline( counter, aone, 90-atwo, cylength, currpos );
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(3);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      zcurr   := basesdev*normaldeviate;
      xcurr   := uniformdeviate(basen*basestep);
      ycurr   := uniformdeviate(basen*basestep);
      aone    := uniformdeviate(180);
      atwo    := 0;
      vcz     := cosd( atwo );
      vcy     := sind( atwo )*sind( aone );
      vcx     := sind( atwo )*cosd( aone );
      currpos := ( xcurr, ycurr, zcurr );
      lenvec  := cylength*( vcx, vcy, vcz ); 
      counter := incr( counter );
      generatedirline( counter, aone, 90-atwo, cylength, currpos );
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

beginfig(2);  
  counter    := 0;
  for i=1 upto basen:
    for j=1 upto basen:
      zcurr   := basesdev*normaldeviate;
      xcurr   := uniformdeviate(basen*basestep);
      ycurr   := uniformdeviate(basen*basestep);
      aone    := uniformdeviate(180);
      atwo    := angsdev*normaldeviate;
      vcz     := cosd( atwo );
      vcy     := sind( atwo )*sind( aone );
      vcx     := sind( atwo )*cosd( aone );
      currpos := ( xcurr, ycurr, zcurr );
      lenvec  := cylength*( vcx, vcy, vcz ); 
      counter := incr( counter );
      generatedirline( counter, aone, 90-atwo, cylength, currpos );
    endfor;
  endfor;
  NL := counter;
  director_invisible( true, thickn, false );
endfig;

end;
