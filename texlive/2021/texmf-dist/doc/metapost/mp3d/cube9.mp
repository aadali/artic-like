% Cube projections
% D. Roegel 23 January 2001
% roegel@loria.fr


input 3danim
output_res:=72;
drawing_scale:=2cm;

projection_type:=1;

%set_plane_(projection_plane)(point_null,vec_I,vec_J);
%set_plane_(projection_plane)(point_null,vec_I,vec_K);

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
  draw_lines(1,2,4,3,1);
  draw_line(2,6);draw_line(4,8);draw_line(1,5);
  draw_lines(5,6,8);
  drawoptions(dashed evenly);
  draw_line(3,7);draw_lines(8,7,5);
  drawoptions();
enddef;

assign_obj("cube","C");

% Trimetric projection
i:=67;
beginfig(167);
  % Observator
  set_point_(Obs)(20*cosd(3.6*i),20*sind(3.6*i),6);
  Obs_phi:=90;Obs_dist:=2;point_of_view_obj("cube",8,Obs_phi);
  draw_obj("cube"); 
endfig;

% Isometric projection
isometric_projection(vec_I,vec_J,vec_K,pnt_obj("cube",8),Obs_dist,90);
  
beginfig(168);
  draw_obj("cube"); 
endfig;

% Dimetric projection  
dimetric_projection(vec_I,vec_J,vec_K,0.5,pnt_obj("cube",8),3*Obs_dist,90);
  
beginfig(169);
  draw_obj("cube"); 
endfig;
 
end.





