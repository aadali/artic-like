% D. Roegel January 2001
% roegel@loria.fr
% Oblique projections of the cube

input 3danim
output_res:=72;
drawing_scale:=2cm;

%set_plane_(projection_plane)(point_null,vec_I,vec_J);
set_plane_(projection_plane)(point_null,vec_I,vec_K);

new_plane_(horizontal_plane)(point_null,vec_I,vec_J);

vardef def_C(expr inst)=
  new_obj_points(inst,8); % 8 points
  set_C_points(inst);
enddef;

vardef set_C_points(expr inst)=
  set_point(1)(0,0,0);set_point(2)(0,0,1);
  set_point(3)(0,1,0);set_point(4)(0,1,1);
  set_point(5)(1,0,0);set_point(6)(1,0,1);
  set_point(7)(1,1,0);set_point(8)(1,1,1);
enddef;

vardef draw_C(expr inst)=
  draw_line(3,7);draw_lines(8,7,5);
  draw_line(2,6);draw_line(4,8);
  draw_lines(5,6,8);draw_lines(4,3);draw_lines(2,4);
  drawoptions(dashed evenly);
  draw_lines(1,2);draw_lines(3,1);draw_line(1,5);
  drawoptions();
enddef;

assign_obj("cube","C");
point_of_view_obj("cube",8,Obs_phi);

oblique_projection(horizontal_plane)(pnt_obj("cube",8),CAVALIER,30,50); % 30

beginfig(100);
  Obs_dist:=2;
  draw_obj("cube"); 
endfig;
  
oblique_projection(horizontal_plane)(pnt_obj("cube",8),CABINET,30,20); % 30
  
beginfig(101);
  Obs_dist:=2;
  draw_obj("cube"); 
endfig;

oblique_projection(horizontal_plane)(pnt_obj("cube",8),30,20,20); % 20

beginfig(102);
  Obs_dist:=2;
  draw_obj("cube"); 
endfig;

new_point(hp);
set_point_(hp)(0,1,1);

set_plane_(horizontal_plane)(hp,vec_I,vec_K);
% such a plane and a second angle of 90 degrees produces a
% ``planometric'' view.
oblique_projection(horizontal_plane)(pnt_obj("cube",8),30,90,20); % 20

vardef draw_C(expr inst)=
  draw_lines(8,7);draw_line(3,7);draw_lines(7,5); 
  draw_lines(1,2);draw_lines(3,1);
  draw_line(1,5);
  draw_lines(4,3);draw_lines(2,4);draw_line(4,8);
  drawoptions(dashed evenly);
  draw_lines(5,6);
  draw_line(2,6);draw_lines(6,8);
  drawoptions();
enddef;

beginfig(103);
  Obs_dist:=2;
  draw_obj("cube"); 
endfig;

end.
  




