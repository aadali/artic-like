% D. Roegel, 7 January 2001
% roegel@loria.fr

verbatimtex
%&latex
\documentclass{article}
\usepackage{mathpple}
\begin{document}
etex

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=10cm;

numeric lsa,lsb,lsc,aasc,aasb,absc;
lsa=9;lsb=8;lsc=4;aasc=60;aasb=40;absc=30;

vardef def_T(expr inst)=
  new_obj_points(inst,19); % 19 points
  set_T_points(inst);
enddef;

vardef set_T_points(expr inst)=
  save hp,kp,sa,sc,inter;
  set_point(1)(0,0,0); % S
  set_point(2)(lsa,0,0); % A
  set_point(4)(lsc*cosd(aasc),lsc*sind(aasc),0); % C
  vec_diff(5,2,1);
  vec_unit(5,5);vec_mult(5,5,lsb*cosd(aasb));
  vec_sum(5,5,1);
  vec_diff(6,4,1);
  vec_unit(6,6);vec_mult(6,6,lsb*cosd(absc));
  vec_sum(6,6,1);
  new_plane(hp)(1,1,1);
  new_plane(kp)(1,1,1);
  new_line(sa)(1,2);
  new_line(sc)(1,4);
  new_line(inter)(1,1);
  def_orth_pl_l_p_(hp)(sa)(pnt(5));
  def_orth_pl_l_p_(kp)(sc)(pnt(6));
  if def_inter_l_pl_pl(inter)(hp)(kp): % there is an intersection
    if def_point_at(3,-lsb,1)(inter):
    else:
      message "PROBLEM";
    fi;
    new_plane(sac)(1,2,4);
    if def_inter_p_l_pl(13)(inter)(sac): else: message "PROBLEM";fi;
    free_plane(sac);
  else:
    message "PROBLEM";
    set_point(3)(1,1,1);
  fi;
  def_right_angle(7,8,9,5,1,3);
  def_right_angle(10,11,12,6,1,3);
  def_right_angle(14,15,16,13,6,3);
  def_right_angle(17,18,19,13,5,3);
  free_line(inter);free_line(sc);free_line(sa);
  free_plane(kp);free_plane(hp);
enddef;

vardef draw_T(expr inst)=
  draw_line(1,2);draw_line_extra(1,4)(-.1,2);draw_line(1,3);
  draw_line(2,3);draw_line(3,4);
  draw_line_extra(3,13)(-0.5,1);
  drawoptions(dashed evenly);
  draw_line(2,4);draw_line(3,5);draw_line(3,6);
  draw_line_extra(3,13)(1,2);
  draw_line(13,6);draw_line(13,5);
  drawoptions();
  draw_double_right_angle(7,8,9,5);
  draw_double_right_angle(10,11,12,6);
  draw_double_right_angle(14,15,16,13);
  draw_double_right_angle(17,18,19,13);

  label_obj.bot(btex $S$ etex,1);
  label_obj.rt(btex $A$ etex,2);
  label_obj.ulft(btex $B$ etex,3);
  label_obj.llft(btex $C$ etex,4);
  
enddef;

square_angle_size:=0.35;
assign_obj("tetra","T");


%for i:=0 upto 20:
  i=0;
  beginfig(100+i);
    % Observator
    set_point_(Obs)(20*cosd(3.6*i-150),20*sind(3.6*i-150),6);
    Obs_phi:=90;Obs_dist:=2;point_of_view_obj("tetra",1,Obs_phi);
    draw_obj("tetra"); 
  endfig;
%endfor;

end.




