% loglogpaper.mp
% L. Nobre G.
% 2012

prologues := 3;

vardef logten( expr Arg )=
  (1.696463*mlog( Arg )/1000)
enddef;

message "log(100) = " & decimal( logten( 100.00 ) );
message "log(300) = " & decimal( logten( 300.00 ) );
message "log(0.01)= " & decimal( logten(   0.01 ) );
message "log(0.03)= " & decimal( logten(   0.03 ) );

beginfig(1);
  numeric xdecs, ydecs, xdecsiz, ydecsiz, llcornx, llcorny, i, j;
  numeric tmpvalx, tmpvaly;
  color colormain, coloraux;
  picture singledec;
  pen penmain, penaux;
  xdecs = 7;
  ydecs = 5;
  xdecsiz = 22mm;
  ydecsiz = 22mm;
  llcornx = 20mm;
  llcorny = 10mm;
  colormain = 0.05white;
  coloraux = blue+green;
  penmain = pencircle scaled 0.26mm;
  penaux = pencircle scaled 0.111mm;
  drawoptions( withpen penaux withcolor coloraux );
  for i=0 upto 10:
    tmpvalx := xdecsiz*(logten(10+i*2)-1);
    draw (tmpvalx,0)--(tmpvalx,ydecsiz);
    tmpvaly := ydecsiz*(logten(10+i*2)-1);
    draw (0,tmpvaly)--(xdecsiz,tmpvaly);
  endfor;
  drawoptions( withpen penmain withcolor colormain );
  for i=1 upto 10:
    tmpvalx := xdecsiz*(logten(10*i)-1);
    draw (tmpvalx,0)--(tmpvalx,ydecsiz);
    tmpvaly := ydecsiz*(logten(10*i)-1);
    draw (0,tmpvaly)--(xdecsiz,tmpvaly);
  endfor;
  singledec = currentpicture;
  currentpicture := nullpicture;
  drawoptions();
  for i=1 upto xdecs:
    for j=1 upto ydecs:
      draw singledec shifted (llcornx+(i-1)*xdecsiz,llcorny+(j-1)*ydecsiz);
    endfor;
  endfor;
  label.bot("0.1",(llcornx,llcorny));
  label.bot("1",(llcornx+xdecsiz,llcorny));
  label.bot("10",(llcornx+2*xdecsiz,llcorny));
  label.bot("100",(llcornx+3*xdecsiz,llcorny));
  label.bot("1000",(llcornx+4*xdecsiz,llcorny));
  label.bot("10000",(llcornx+5*xdecsiz,llcorny));
  label.bot("100000",(llcornx+6*xdecsiz,llcorny));
  label.bot("1000000",(llcornx+7*xdecsiz,llcorny));
  label.lft("0.0001",(llcornx,llcorny));
  label.lft("0.001",(llcornx,llcorny+ydecsiz));
  label.lft("0.01",(llcornx,llcorny+2*ydecsiz));
  label.lft("0.1",(llcornx,llcorny+3*ydecsiz));
  label.lft("1",(llcornx,llcorny+4*ydecsiz));
  label.lft("10",(llcornx,llcorny+5*ydecsiz));
endfig;

end.
  