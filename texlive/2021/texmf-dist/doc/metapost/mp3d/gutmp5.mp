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
  new_obj_points(inst,4); % 4 points
  set_T_points(inst);
enddef;

vardef set_T_points(expr inst)=
  save hp,kp,sa,sc,inter;new_point(h);new_point(k);
  set_point(1)(0,0,0); % S
  set_point(2)(lsa,0,0); % A
  set_point(4)(lsc*cosd(aasc),lsc*sind(aasc),0); % C
  vec_diff_(h,pnt(2),pnt(1));
  vec_unit_(h,h);vec_mult_(h,h,lsb*cosd(aasb));
  vec_sum_(h,h,pnt(1));
  vec_diff_(k,pnt(4),pnt(1));
  vec_unit_(k,k);vec_mult_(k,k,lsb*cosd(absc));
  vec_sum_(k,k,pnt(1));
  new_plane(hp)(1,1,1);
  new_plane(kp)(1,1,1);
  new_line(sa)(1,2);
  new_line(sc)(1,4);
  new_line(inter)(1,1);
  def_orth_pl_l_p_(hp)(sa)(h);
  def_orth_pl_l_p_(kp)(sc)(k);
  if def_inter_l_pl_pl(inter)(hp)(kp): % there is an intersection
    if def_point_at(3,-lsb,1)(inter):
    else:
      message "PROBLEM";
    fi;
  else:
    message "PROBLEM";
    set_point(3)(1,1,1);
  fi;
  free_line(inter);free_line(sc);free_line(sa);
  free_plane(kp);free_plane(hp);
  free_point(k);free_point(h);
enddef;

vardef draw_T(expr inst)=
  draw_line(1,2);draw_line(1,4);draw_line(1,3);
  draw_line(2,3);draw_line(3,4);
  drawoptions(dashed evenly);
  draw_line(2,4);
  drawoptions();
  label_obj.bot(btex $S$ etex,1);
  label_obj.rt(btex $A$ etex,2);
  label_obj.top(btex $B$ etex,3);
  label_obj.lft(btex $C$ etex,4);
  
enddef;

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




