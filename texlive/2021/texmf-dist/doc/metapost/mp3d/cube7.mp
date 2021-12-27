% The nine vanishing points of a cube
% D. Roegel, 22 January 2001
% Warning: this figure is very sensitive to small motions!

input 3danim
input 3dgeom
output_res:=72;
drawing_scale:=1cm;

def draw_arrowline_extra(expr i,j)(expr exi,exj)=
  drawarrow exi[z[ipnt_(i)],z[ipnt_(j)]]--exj[z[ipnt_(i)],z[ipnt_(j)]];
enddef;


vardef def_C(expr inst)=
  new_obj_points(inst,17); % 17 points
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
  set_line(l)(5,3);
  if not def_vanishing_point_p_l_pl(13)(l)(screen):
    message "no vanishing point 5";
    set_point(13)(0,0,0);
  fi;
  set_line(l)(5,8);
  if not def_vanishing_point_p_l_pl(14)(l)(screen):
    message "no vanishing point 6";
    set_point(14)(0,0,0);
  fi;
  set_line(l)(6,7);
  if not def_vanishing_point_p_l_pl(15)(l)(screen):
    message "no vanishing point 7";
    set_point(15)(0,0,0);
  fi;
  set_line(l)(5,2);
  if not def_vanishing_point_p_l_pl(16)(l)(screen):
    message "no vanishing point 8";
    set_point(16)(0,0,0);
  fi;
  set_line(l)(6,1);
  if not def_vanishing_point_p_l_pl(17)(l)(screen):
    message "no vanishing point 9";
    set_point(17)(0,0,0);
  fi;
  
  free_line(l);
enddef;

vardef draw_C(expr inst)=
  draw_lines(1,2,4,3,1);draw_lines(5,6,8,7,5);
  draw_line(1,5);draw_line(2,6);draw_line(3,7);draw_line(4,8);
  pickup pencircle scaled 4pt;
  drawoptions(withcolor red);
  for i:=9 upto 17:
    if z[ipnt_(i)]=(10,10):draw_line(i,i);
    fi;
  endfor;
  drawoptions();
  pickup pencircle scaled 1pt;
  for i:=9 upto 17:draw_line(i,i);endfor;
  pickup pencircle scaled 0.4pt;
  drawoptions(dashed evenly scaled 0.5);
  for j:=1,3,5,7:draw_line(j,9);endfor;
  for j:=3,4,7,8:draw_line(j,10);endfor;
  for j:=1 upto 4:draw_line(j,11);endfor;
  draw_line(1,12);draw_line(2,12);
  drawoptions(dashed evenly scaled 0.5 withcolor red);
  for j:=5,6:draw_line(j,13);endfor;
  drawoptions(dashed evenly scaled 0.5 withcolor green);
  for j:=1,5:draw_line(j,14);endfor;
  drawoptions(dashed evenly scaled 0.5 withcolor blue);
  for j:=2,6:draw_line(j,15);endfor;
  drawoptions(dashed evenly scaled 0.5 withcolor green);
  for j:=2,4:draw_line(j,16);endfor;
  drawoptions(dashed evenly scaled 0.5 withcolor red);
  for j:=6,8:draw_line(j,17);endfor;

  drawoptions();
enddef;

new_plane_(screen)(0,0,0); % must be defined before |assign_obj|
                           % because |assign_obj| calls |set_C_point|

% we change the field, because we are very close
h_field:=170; % degrees 
v_field:=170; % degrees


assign_obj("cube","C");

%for i:=0 upto 100:
i=96;
  beginfig(100+i);
    % Observator
    set_point_(Obs)(1.7*cosd(3.6*i),1.8*sind(3.6*i),1.8);
    Obs_phi:=90;Obs_dist:=0.5;point_of_view_obj("cube",8,Obs_phi);
    reset_obj("cube");
    draw_obj("cube");
  endfig;
%endfor; 

end.




