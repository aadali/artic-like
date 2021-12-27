% D. Roegel
% Three orthogonal projections of a tetrahedron
% January 7, 2001
% roegel@loria.fr

verbatimtex
%&latex
\documentclass{article}
\usepackage{mathpple}
\begin{document}
etex



input 3dgeom;   % definitions for geometry
input 3danim;   % animations

% We modify some of the default settings.
output_res:=72;
drawing_scale:=7cm;
filled_faces:=false;
show_animation_parameters:=true;

%projection_type:=1; % parallel projection

%================ DEFINITION OF OBJECT (BEGIN) =========================
vardef def_fig(expr inst)=
  new_obj_points(inst,28); % 28 points 4+3*(4+4)
  set_fig_points(inst);
enddef;

% Definition of the points that belong to the general figure.
vardef set_fig_points(expr inst)=
  set_point(1)(0.3,0.45,0.8); % A
  set_point(2)(0.6,0.45,0.8); % B
  set_point(3)(0.5,0.6,0.5);  % C
  set_point(4)(0.5,0.3,0.5);  % D

  % the three planes (four corners at a time)
  % horizontal
  set_point(5)(0,0,-.2);
  set_point(6)(1,0,-.2);
  set_point(7)(1,1,-.2);
  set_point(8)(0,1,-.2);
  % vertical left
  set_point(9) (-.45,0,0);
  set_point(10)(-.45,1,0);
  set_point(11)(-.45,1,1);
  set_point(12)(-.45,0,1);
  % vertical front
  set_point(13)(0,-.45,0);
  set_point(14)(1,-.45,0);
  set_point(15)(1,-.45,1);
  set_point(16)(0,-.45,1);

  % we now compute the verticals and therefore introduce plane structures
  new_plane(f)(9,10,11);
  new_plane(g)(13,14,15);
  new_plane(h)(5,6,7);

  % bottow figure
  def_vert_pl(17)(1)(h);
  def_vert_pl(18)(2)(h);
  def_vert_pl(19)(3)(h);
  def_vert_pl(20)(4)(h);

  % left figure
  def_vert_pl(21)(1)(f);
  def_vert_pl(22)(2)(f);
  def_vert_pl(23)(3)(f);
  def_vert_pl(24)(4)(f);

  % front figure
  def_vert_pl(25)(1)(g);
  def_vert_pl(26)(2)(g);
  def_vert_pl(27)(3)(g);
  def_vert_pl(28)(4)(g);

  free_plane(h);free_plane(g);free_plane(f);
enddef;

% Method to draw a general figure
vardef draw_fig(expr inst)=
  % tetrahedron
  pickup pencircle scaled 1pt;
  drawoptions();
  draw_line(1,2);draw_line(1,3);draw_line(2,3);
  draw_line(2,4);draw_line(3,4);
  drawoptions(dashed withdots);draw_line(1,4);
  drawoptions();

  % projection on f:
  draw_line(21,22);draw_line(21,23);draw_line(21,24);
  draw_line(22,23);draw_line(22,24);draw_line(23,24);

  % projection on g:
  draw_line(25,26);draw_line(25,27);draw_line(25,28);
  draw_line(26,27);draw_line(26,28);draw_line(27,28);

  % projection on h:
  draw_line(17,18);draw_line(17,19);draw_line(17,20);
  draw_line(18,19);draw_line(18,20);draw_line(19,20);

  pickup pencircle scaled .4pt;

  % les plans
  draw_line(5,6);draw_line(6,7);draw_line(7,8);draw_line(8,5);
  draw_line(9,10);draw_line(10,11);draw_line(11,12);draw_line(12,9);
  draw_line(13,14);draw_line(14,15);draw_line(15,16);draw_line(16,13);
  pickup pencircle scaled .4pt;

  % dotted projection lines:
  save i;drawoptions(dashed evenly);
  for i:=1 upto 4: draw_line(i,i+16);endfor;
  for i:=1 upto 4: draw_line(i,i+20);endfor;
  for i:=1 upto 4: draw_line(i,i+24);endfor;

  % and a few labels:
  label_obj.top(btex $A$ etex,1);
  label_obj.llft(btex $B$ etex,2);
  label_obj.lrt(btex $C$ etex,3);
  label_obj.llft(btex $D$ etex,4);

  label_obj.bot(btex $H$ etex,7);
  label_obj.top(btex $F$ etex,12);
  label_obj.llft(btex $G$ etex,14);

enddef;

assign_obj("figa","fig");

%================ DEFINITION OF OBJECT (END) =========================

% Call of drawing instruction
%for i:=0 upto 20:%99:
i=12;
  beginfig(100+i);
    % Positionning of the observer
    set_point_(Obs)(4*cosd(3.6*i),4*sind(3.6*i),2.2);
    Obs_phi:=90;Obs_dist:=2;point_of_view_obj("figa",1,Obs_phi);
    draw_obj("figa"); % main figure
  endfig;
%endfor; show_animation_bbox;
end.





