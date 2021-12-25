% D. Roegel, 9 January 2001
% roegel@loria.fr

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=15cm;

vardef def_T(expr inst)=
  new_obj_points(inst,12); % 12 points
  set_T_points(inst);
enddef;

vardef set_T_points(expr inst)=
  set_point(1)(0,0,0);
  set_point(2)(1,0,0);
  set_point(3)(1,1,0);
  set_point(4)(0,1,0);
  set_point(5)(0,0,1);
  set_point(6)(1,0,1);
  set_point(7)(1,1,1);
  set_point(8)(0,1,1);
  new_vec(v);
  vec_diff_(v,pnt(1),pnt(4));
  vec_sum_(pnt(9),pnt(1),v);
  mid_point(10,1,2);
  vec_sum_(pnt(11),pnt(2),v);
  mid_point(12,7,8);
enddef;

vardef draw_T(expr inst)=
  draw_line(1,2);draw_line(2,3);
  draw_line(5,6);draw_line(6,7);draw_line(7,8);draw_line(8,5);
  draw_line(1,5);draw_line(2,6);draw_line(3,7);
  draw_line(9,10);draw_line(10,11);draw_line(11,9);
  draw_line(5,12);draw_line(12,6);draw_line(6,5);
  drawoptions(dashed evenly);
  draw_line(3,4);draw_line(4,1);draw_line(4,8);
  draw_line(5,9);draw_line(10,12);draw_line(6,11);
  drawoptions();
  fill z[ipnt_(9)]--z[ipnt_(10)]--z[ipnt_(11)]--cycle withcolor 0.8white;
  fill z[ipnt_(5)]--z[ipnt_(12)]--z[ipnt_(6)]--cycle withcolor 0.8white;
  draw_equal_marks(1,4,2);
  draw_equal_marks(1,5,2);
enddef;

assign_obj("tetra","T");

%for i:=0 upto 30:
i=23;
  beginfig(100+i);
    % Observator
    set_point_(Obs)(10*cosd(3.6*i-150),10*sind(3.6*i-150),3);
    Obs_phi:=90;Obs_dist:=2;point_of_view_obj("tetra",1,Obs_phi);
    draw_obj("tetra"); 
  endfig;
%endfor;

end.




