% D. Roegel January 2001
% roegel@loria.fr

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=2cm;

vardef def_C(expr inst)=
  new_obj_points(inst,12); % 12 points
  set_C_points(inst);
enddef;

vardef set_C_points(expr inst)=
  set_point(1)(0,0,0);set_point(2)(0,0,1);
  set_point(3)(0,1,0);set_point(4)(0,1,1);
  set_point(5)(1,0,0);set_point(6)(1,0,1);
  set_point(7)(1,1,0);set_point(8)(1,1,1);
  % vanishing points
  def_screen_pl(screen);
  new_line(l)(1,2);
  if not def_vanishing_point_p_l_pl(9)(l)(screen):
    message "no vanishing point 1";
    set_point(9)(0,0,0);
  fi;
  set_line(l)(1,3);
  if not def_vanishing_point_p_l_pl(10)(l)(screen):
    message "no vanishing point 2";
    set_point(10)(0,0,0);
  fi;
  set_line(l)(1,5);
  if not def_vanishing_point_p_l_pl(11)(l)(screen):
    message "no vanishing point 3";
    set_point(11)(0,0,0);
  fi;
  set_line(l)(1,7);
  if not def_vanishing_point_p_l_pl(12)(l)(screen):
    message "no vanishing point 4";
    set_point(12)(0,0,0);
  fi;
  free_line(l);
enddef;

vardef draw_C(expr inst)=
  draw_lines(1,2,4,3,1);draw_lines(5,6,8,7,5);
  draw_line(1,5);draw_line(2,6);draw_line(3,7);draw_line(4,8);
  pickup pencircle scaled 4pt;
  drawoptions(withcolor red);
  if z[ipnt_(9)]=(10,10):draw_line(9,9);
  fi;
  if z[ipnt_(10)]=(10,10):draw_line(10,10);
  fi;
  if z[ipnt_(11)]=(10,10):draw_line(11,11);
  fi;
  if z[ipnt_(12)]=(10,10):draw_line(12,12);
  fi;
  drawoptions();
  pickup pencircle scaled 1pt;
  draw_line(9,9);draw_line(10,10);draw_line(11,11);draw_line(12,12);
  pickup pencircle scaled 0.4pt;
  drawoptions(dashed evenly scaled 0.5);
  for j:=1,3,5,7:draw_line(j,9);endfor;
  for j:=3,4,7,8:draw_line(j,10);endfor;
  for j:=1 upto 4:draw_line(j,11);endfor;
  draw_line(1,12);draw_line(2,12);
  drawoptions();
enddef;

new_plane_(screen)(0,0,0); % must be defined before |assign_obj|
                           % because |assign_obj| calls |set_C_point|

% on change le champ car on est très près
h_field:=160; % degrees 
v_field:=160; % degrees


assign_obj("cube","C");

%for i:=0 upto 100:
i=97;
  beginfig(100+i);
    % Observator
    set_point_(Obs)(1.8*cosd(3.6*i),1.8*sind(3.6*i),2);
    Obs_phi:=90;Obs_dist:=0.5;point_of_view_obj("cube",8,Obs_phi);
    reset_obj("cube");
    draw_obj("cube");
  endfig;
%endfor; 

end.




