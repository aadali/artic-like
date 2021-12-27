verbatimtex
%&latex
\documentclass{article}

\begin{document}
etex

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=20cm;

vardef def_C(expr inst)=
  new_obj_points(inst,21); % 21 points
  set_C_points(inst);
enddef;

vardef set_C_points(expr inst)=
  set_point(1)(0,0,0);set_point(2)(0,0,1);
  set_point(3)(0,1,0);set_point(4)(0,1,1);
  set_point(5)(1,0,0);set_point(6)(1,0,1);
  set_point(7)(1,1,0);set_point(8)(1,1,1);

  % corners of the projection plane
  %set_point(9)(12.92993,12.37372,6.83334);
  %set_point(10)(15.68588,8.48413,6.8334);
  %set_point(11)(16.27899,8.90436,4.12851);
  %set_point(12)(13.52304,12.79395,4.12845);

  %set_point(9)(13.76718,11.50638,6.15714);
  %set_point(10)(15.14517,9.56158,6.15717);
  %set_point(11)(15.44174,9.7717,4.80472);
  %set_point(12)(14.06375,11.71649,4.80469);

  set_point(9)(0.80682,2.27571,1.939);
  set_point(10)(2.03946,0.2357,1.93903);
  set_point(11)(2.46796,0.49461,0.63113);
  set_point(12)(1.23532,2.53462,0.6311);

  % shown observator
  set_point(13)(3.23608,2.35114,2);

  % intersections with the projection plane
  new_plane(pl)(9,10,11);
  new_line(l)(1,13);
  if not def_inter_p_l_pl(14)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(2,13);
  if not def_inter_p_l_pl(15)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(3,13);
  if not def_inter_p_l_pl(16)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(4,13);
  if not def_inter_p_l_pl(17)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(5,13);
  if not def_inter_p_l_pl(18)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(6,13);
  if not def_inter_p_l_pl(19)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(7,13);
  if not def_inter_p_l_pl(20)(l)(pl):
    message "PROBLEM ";
  fi;
  set_line(l)(8,13);
  if not def_inter_p_l_pl(21)(l)(pl):
    message "PROBLEM ";
  fi;
  
  free_line(l);
  free_plane(pl);
  
enddef;

vardef draw_C(expr inst)=
  draw_lines(1,2,4,3,1);draw_lines(5,6,8,7,5);
  draw_line(1,5);draw_line(2,6);draw_line(3,7);draw_line(4,8);
  draw_lines(9,10,11,12,9);
  pickup pencircle scaled 2pt;
  draw_line(13,13);draw_line(1,1);draw_line(2,2);draw_line(3,3);draw_line(7,7);
  pickup pencircle scaled .4pt;
  draw fullcircle scaled 6pt shifted z[ipnt_(8)];
  label_obj.lft(btex \textit{Obs}$(3.23,2.35,2)$ etex,13);
  label_obj.rt(btex $(0,0,0)$ etex,1);
  label_obj.rt(btex $(0,0,1)$ etex,2);
  label_obj.bot(btex $(0,1,0)$ etex,3);
  label_obj.llft(btex $(1,1,0)$ etex,7);
  
  drawoptions(dashed withdots);
  draw_line(1,13);draw_line(2,13);draw_line(3,13);draw_line(4,13);
  draw_line(5,13);draw_line(6,13);draw_line(7,13);draw_line(8,13);
  drawoptions();
  draw_line(14,15,17,16,14);draw_lines(18,19,21,20,18);
  draw_line(14,18);draw_line(15,19);draw_line(16,20);draw_line(17,21);
enddef;

assign_obj("cube","C");

for i:=28 upto 28:
  beginfig(100+i);
    % Observator
    set_point_(Obs)(20*cosd(3.6*i),20*sind(3.6*i),4);
    show_point("Obs",Obs);
    Obs_phi:=90;Obs_dist:=2;point_of_view_obj("cube",14,Obs_phi);
    draw_obj("cube"); 
  endfig;
endfor; 

end.







