% mysquaresectionbar.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

% bash: export TEX=latex
% tcsh: setenv TEX latex

prologues := 1;

f := (3.5,3,2.5);
Spread := 175;
beginfig(1);
  numeric alfa, beta;
  alfa = 22.62;
  kindofcube(false,false,(-0.65,0,-0.6),  90,0,0,1.35,0.05,0.05);
  kindofcube(false,false,(-0.65,0,0.6),   90,0,0,0.35,0.05,0.05);
  kindofcube(false,false,(-0.6,-0.4,0),   90,90-alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(-0.6,0.4,0),    90,90+alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(-0.6,-0.05,0.05),0,90,0,1.35,0.05,0.05);
  kindofcube(false,false,(0.6,-0.05,0.05), 0,90,0,1.35,0.05,0.05);
  kindofcube(false,false,(0,0,-0.3),       0,alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(0,0,0.3),        0,-alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(0.6,-0.4,0),    90,90-alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(0.6,0.4,0),     90,90+alfa,0,1.35,0.05,0.05);
  kindofcube(false,false,(0.65,0,-0.6),   90,0,0,1.35,0.05,0.05);
  kindofcube(false,false,(0.65,0,0.6),    90,0,0,0.35,0.05,0.05);
  kindofcube(false,false,(0,0,0.65),       0,0,0,1.35,0.05,0.05);
endfig;

end;
  