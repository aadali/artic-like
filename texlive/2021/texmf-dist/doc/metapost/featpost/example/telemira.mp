pair thecenter;
thecenter = (105mm,150mm);

vardef instructio( expr rayc ) =
  draw fullcircle scaled rayc shifted thecenter
enddef;

beginfig(1);
  numeric i, u, maxray;
  u = 1mm;
  maxray = 85u;
  for i=30u step 1u until 40u:
    instructio( i );
  endfor;
  for i=153u step 1u until 165u:
    instructio( i );
  endfor;
  for i=0 step 90 until 270:
    draw thecenter--(thecenter+maxray*dir(i));
  endfor;
endfig;
end.