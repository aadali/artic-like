% bughunt.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

%SphericalDistortion := true;
%f 		    := 0.3*(5,4,2);
%Spread		    := 40;
%ShiftV		    := (200,200);

SphericalDistortion := false;
f 		    := 3.3*(5,4,2);
Spread		    := 70;
ShiftV		    := (250,250);

beginfig(1);
    color nalfa, gammacnt, vertex;
    numeric altura, grids, radius, aux, i, numc, slope;
    
    altura 		= 1;
    radius		= 0.3;
    grids		= 1;
    numc		= 8;
    slope		= 15;

    for i=1 upto numc:
        aux		:= (i-1-0.5*(numc-1))*grids;
        gammacnt 	:= (0,aux,0);
        vertex		:= gammacnt+0.1*(f-gammacnt)+(0,0,aux/slope);
	nalfa           := vertex - gammacnt;
        rigorousdisc( 0, true, gammacnt, radius, nalfa );
    endfor;
endfig;

beginfig(2);
  path cpath;
  color vert;
  numeric i;
  f := (4,0,0.5);
  vert = 0.5*(0,1,1);
  cpath = rigorouscone(true, black, blue, 1, vert);
%   draw rp(vert);
  draw cpath; show length cpath;
  for i=0 upto length cpath:
    dotlabel.top( decimal(i),point i of cpath);
  endfor;
%   cpath := rigorouscircle( black, blue, 1);
%   draw cpath;
%   for i=0 upto length cpath:
%     dotlabel.bot( decimal(i),point i of cpath);
%   endfor;
endfig;

end;


