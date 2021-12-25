% goldenellipse.mp
% L. Nobre G.
% 2012

beginfig(1);
  numeric u, phi, i, dstep, goldang;
  pen thepen;
  path ellia, ellib, ellic, ellid, ellip;
  path evola, evolb, evolc, evold, evolu;
  path squar;
  pair paira, pairb, pairc, paird;
  pair ellfa, squab, ellfc, squad;
  pair goldp;
  u = 5cm;
  phi = 0.5*(1+sqrt(5));
  dstep = 5;
  goldang = angle(1,1/phi);
  thepen = pencircle scaled 0.05u;
  paira = (phi,0);
  pairb = (0,1);
  pairc = (-phi,0);
  paird = (0,-1);
  ellfa = (1,0);
  squab = (0,-phi);
  ellfc = (-1,0);
  squad = (0,phi);
  goldp = (phi*cosd(goldang),sind(goldang));
  ellia = (phi*cosd(dstep),sind(dstep))
    for i=2dstep step dstep until 90-dstep:
    ..(phi*cosd(i),sind(i))
  endfor;
  ellib = reverse ellia xscaled -1;
  ellic = ellia scaled -1;
  ellid = reverse ellia yscaled -1;
  ellip = paira--ellia--pairb--ellib--pairc--ellic--paird--ellid--cycle;
  evola = (cosd(dstep)**3,-phi*(sind(dstep)**3))
    for i=2dstep step dstep until 90-dstep:
    ..(cosd(i)**3,-phi*(sind(i)**3))
  endfor;
  evolb = reverse evola xscaled -1;
  evolc = evola scaled -1;
  evold = reverse evola yscaled -1;
  evolu = ellfa--evola--squab--evolb--ellfc--evolc--squad--evold--cycle;
  squar = paira--squab--pairc--squad--cycle;
  z1 = goldp+whatever*(1,-1);
  y1 = 0;
  fill (x1,0)--(0,x1)--(-x1,0)--(0,-x1)--cycle withcolor blue;
  fill squar withcolor red;
  fill ellip;
  fill ellfa--pairb--ellfc--paird--cycle withcolor red;
  draw evolu withcolor blue;
  picture tmp;
  tmp = currentpicture;
  currentpicture := nullpicture;
  draw tmp rotated 45 scaled u shifted (10.5cm,14cm) withpen thepen;
endfig;
end.