%jose samuel perez huerta
%mayo 2007
%paladium lattice in hydrogen

input featpost3Dplus2D;

% Under the adsorption process,
% the hydrogen molecules are retained in the surface of the
% metal forming a gaseous layer. Nevertheless, the hydrogen molecule is
% too great to get up in the structure of a solid. When the amount of
% hydrogen is increased, the palladium dissociates molecules of the gas
% so that now in its atomic form, they are diffused between the
% interstitial sites of lattice formation of palladium hydrides In this
% process, the lattice undergoes to isotropic expansion, lattice
% constant is increased until a 3\% in each space direction, increasing
% the total volume until in a 10\%

beginfig(1);
  numeric sph, sphc, ind, indos, suma, intres, i, j, coos, size, num;
  color xdelta, ydelta, zdelta, W[], fro, pro;
  path p[], pp[], ppp, pppp;
  f := 1.5*(10,5,3);
%  ParallelProj := true;
  Spread := 50;
  sph := 0.15;
  sphc := 0.07;
  size = 6;
  num=200;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xdelta := (sph,0,0);
  ydelta := (0,sph,0);
  zdelta := (0,0,sph);
  suma:=0;
  for ind=1 upto 3:
    for indos=1 upto 10:
      for intres=1 upto 10:
	suma:=suma+1;
	V[suma] :=0.5*(indos,intres,ind);
	W[suma] :=0.5*((indos-1,intres-1,ind-1)+0.5*(1,1,1));
      endfor;
    endfor;
  endfor;
  
  for ind=1 upto suma:
    p[ind] := rigorousfearpath( V[ind], sph );
    pp[ind] := rigorousfearpath( W[ind], sphc );
  endfor;
  
  for ind=1 upto suma:
    fill pp[ind] withcolor red;
    draw pp[ind] withcolor black;
    fill p[ind] withcolor .75white;
    draw p[ind] withcolor blue;
  endfor;
%%%%%%%%%%%%%%%%%%%%%%  
  for i=1 upto num:
    for j=1 upto 2:
      coos[j] := 1.7*uniformdeviate(size*.5);
    endfor;
    fro := (coos1,coos2, 4*(.5)+(.3)*normaldeviate);
    pro := fro+(0,sph,0);
    ppp := rigorousfearpath( fro, sphc);
    pppp := rigorousfearpath( pro, sphc);
    fill ppp withcolor red;
    draw ppp  withcolor (.5,.4,0);
    fill pppp withcolor red;
    draw pppp  withcolor black;
  endfor;
%%%%%%%%%%%%%%%%%%%%%%%%
endfig;
end;
