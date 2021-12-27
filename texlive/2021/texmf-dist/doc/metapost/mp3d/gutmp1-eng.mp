% D. Roegel, 7 January 2001
% roegel@loria.fr

verbatimtex
%&latex
\documentclass{article}
\begin{document}
etex

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=7cm;

numeric b,c,d,e;b=70;c=0;d=20;e=50; % parameterized angles

vardef def_T(expr inst)=
  new_obj_points(inst,6); % 6 points
  set_T_points(inst);
enddef;

vardef set_T_points(expr inst)=
  set_point(1)(0,0,0);
  set_point(2)(3*cosd(b),3*sind(b),0);
  set_point(3)(6*cosd(c),6*sind(c),0);
  new_vec(v_a);new_vec(v_b);
  vec_def_vec_(v_a,vec_I);
  vec_rotate_(v_a,vec_K,d);
  vec_prod_(v_b,v_a,vec_K);
  vec_rotate_(v_a,v_b,e);
  vec_mult_(v_a,v_a,4.5);
  vec_sum_(pnt(4),pnt(1),v_a);
  free_vec(v_b);free_vec(v_a);
  % Determination of I and J:
  % I=A + unit(AB)
  vec_diff(5,2,1);
  vec_unit(5,5);
  vec_sum(5,5,1);
  % J=A + 4*unit(AC)
  vec_diff(6,3,1);
  vec_unit(6,6);
  vec_mult(6,6,4);
  vec_sum(6,6,1);
  
enddef;

vardef draw_T(expr inst)=
  drawoptions(dashed evenly);
  draw_line(1,2);draw_line(1,3);draw_line(1,4);
  drawoptions();
  draw_line(2,3);draw_line(2,4);draw_line(3,4);
  pickup pencircle scaled 2pt;
  draw_line(5,5);draw_line(6,6);
  label_obj.ulft(btex $A_1$ etex,1);
  label_obj.rt(btex $B_2$ etex,2);
  label_obj.lft(btex $C_3$ etex,3);
  label_obj.top(btex $D_4$ etex,4);
  label_obj.top(btex $I_5$ etex,5);
  label_obj.top(btex $J_6$ etex,6);

enddef;

assign_obj("tetra","T");

%for i:=0 upto 20:
i=8;
  beginfig(100+i);
    % Observator
    set_point_(Obs)(20*cosd(3.6*i),20*sind(3.6*i),6);
    Obs_phi:=90;Obs_dist:=2;point_of_view_obj("tetra",1,Obs_phi);
    draw_obj("tetra"); 
  endfig;
%endfor;

end  
