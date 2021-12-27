% telheiro.mp
% L. Nobre G.
% 2013

input featpost3Dplus2D;

prologues := 1;

f := (6,1,3);

vardef drawpat( expr Ind )( text Nam ) =
  draw rp(Nam[Ind]) for i=Ind-1 downto 1: --rp(Nam[i]) endfor
    for i=1 upto Ind: --rp((-X(Nam[i]),Y(Nam[i]),Z(Nam[i]))) endfor
enddef;

vardef xz( expr C )= (X(C),Z(C)) enddef;

vardef xy( expr C )= (X(C),Y(C)) enddef;

vardef recter( expr Pic ) = 0.5*((llcorner Pic)+(urcorner Pic)) enddef;

beginfig(1);
  numeric  semimajor, semiminor, vertez, vertey, centery;
  color verte, cente, maj, min, ver[], dia[], hor[], el;
  numeric i, cou, ste, axe, maxang;
  semimajor = 1.4;
  semiminor = 1;
  vertez = 1.1;
  vertey = 0.25;
  centery = 0.35;
  ste = 2;
  axe =1.5;
  verte = (0,-vertey,vertez);
  cente = (0,centery,0);
  maj = (0,semimajor,0);
  min = (semiminor,0,0);
  color cen;
  cen = (0,vertey+centery+semimajor,-vertez);
  draw rp(verte)--rp(verte+cen);
  cou = 1;
  i = 1;
  forever:
    el := cente+planarrotation( maj, min, i );
    ver[cou] = whatever[verte,el];
    Y(ver[cou]) = 0;
    dia[cou] = whatever[verte,el];
    Z(dia[cou]) = Y(dia[cou]);
    hor[cou] = whatever[verte,el];
    Z(hor[cou]) = 0;
    exitif (Z(ver[cou])>=vertez) or (Z(ver[cou])<0) or (Y(hor[cou])<0);
    cou := cou+1;
    i := i+ste;
  endfor;
  maxang = i;
  cartaxes(axe,axe,axe);
  draw rp(black)--rp((-axe,0,0));
  drawpat(cou)(ver);
  drawpat(cou)(hor);
  drawoptions( withcolor blue );
  drawpat(cou)(dia);
endfig;

beginfig(2);
  numeric u;
  u = 5cm;
  draw (xz(ver[cou]) for i=cou-1 downto 1: --xz(ver[i]) endfor
    for i=1 upto cou: --(-X(ver[i]),Z(ver[i])) endfor
    for i=cou downto 1: --(-X(hor[i]),Y(hor[i])) endfor
    for i=1 upto cou: --xy(hor[i]) endfor --cycle) scaled u;
  picture surface;
  surface = currentpicture;
  currentpicture := nullpicture;
  numeric ang, v, d, h;
  pair ve[], di[], ho[], vec, possurf, posbase;
  possurf = recter( surface );
  ang = 0;
  for i=1 upto cou:
    ang:= ang+getangle( cen, hor[i]-verte );
    cen:= hor[i]-verte;
    v := conorm(verte-ver[i]);
    d := conorm(verte-dia[i]);
    h := conorm(cen);
    vec := dir(90-ang);
    ve[i] = u*v*vec;
    di[i] = u*d*vec;
    ho[i] = u*h*vec;    
  endfor;
  draw di1 for i=2 upto cou: --di[i] endfor withcolor blue;
  draw for i=cou downto 1: ve[i]-- endfor for i=1 upto cou: ho[i]-- endfor cycle;
  draw currentpicture xscaled -1;
  posbase = recter( currentpicture );
  draw surface shifted (posbase-possurf);
  produce_auto_scale;
endfig;
  
end.