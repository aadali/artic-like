% brownellips.mp
% L. Nobre G.
% IYP (2005)

input featpost3Dplus2D;

numeric ParameterA, ParameterB, ParameterC;
ParameterA = 0.5;            % f
ParameterB = 0.5;            % viewcentr
ParameterC = 405;            % main arm angle

f := ParameterA[6*(3,4,2),4*(2,-3,1)];
viewcentr := ParameterB[black,(-2,-3,1)];
Spread := 16;

beginfig(0);
  numeric sray, armb, mara, flen, pif, pers, perb, rads, radb, pert, j;
  numeric nums, lens, i, bray, unan, stan, acta, arms, mang, actl, coun;
  numeric rrat, auxd, anga, sqsd, k;
  color c[][], vr, vt, vi, vb, auxc;
  path elo[];
  string inde, comm;
  stan = 0;
  nums = 30;
  sray = 0.95;
  bray = 2;                  % bray must be 2*sray if you want ellipses
  armb = 5;
  arms = 2.5;
  sqsd = 0.3;
  pif = 3.14159/180;
  rrat = 1-bray/sray;
  flen = armb +-+ (bray-sray);
  mara = angle( bray-sray, flen );
  rads = pif*2*mara;
  radb = pif*(360-2*mara);
  pers = rads*sray;
  perb = radb*bray;
  pert = pers + perb + 2*flen;
  lens = pert/nums;
  unan = lens/bray/pif;
  acta = stan + ParameterC;
  mang = acta mod unan;
  drawoptions( withcolor 0.7white ); 
  setthestage( 20, 2*(armb+2*arms) );
  drawoptions();
  draw
  rp((armb*cosd(stan)+arms*cosd(rrat*stan),
      armb*sind(stan)+arms*sind(rrat*stan),
      0))
  for i=stan+2 step 2 until acta:
    ...rp((armb*cosd(i)+arms*cosd(rrat*i),armb*sind(i)+arms*sind(rrat*i),0))
  endfor;
  auxc:=(armb*cosd(acta)+arms*cosd(rrat*acta),
    armb*sind(acta)+arms*sind(rrat*acta),
    0);
  verygoodcone(false,auxc+2*sqsd*blue,-blue,sqsd,auxc);
  vb := (armb*cosd(acta),armb*sind(acta),2*sqsd);
  vi := 0.5*sqsd*(sind(rrat*acta),-cosd(rrat*acta),0);
  kindofcube(false,true,vb+vi,rrat*acta,0,0,arms,sqsd,sqsd);
  coun = 0;
  for i=unan-mang+mara-180 step unan until 180-mara:
    coun := incr( coun );
    vr := (-bray*cosd(i),-bray*sind(i),3*sqsd);
    vt := 0.25*lens*(-sind(i),cosd(i),0);
    c[coun]1 = vt+vr;
    c[coun]2 = vt+vr+lens*blue;
    c[coun]3 = vr-vt+lens*blue;
    c[coun]4 = vr-vt;
    c[coun]5 = 0.25*(c[coun]1+c[coun]2+c[coun]3+c[coun]4);
    j := i;
  endfor;
  auxd = lens-pif*(180-mara-j)*bray;
  vi := (-bray*cosd(180-mara),-bray*sind(180-mara),3*sqsd);
  for i=auxd step lens until flen:
    coun := incr( coun );
    vr := vi + (i*sind(mara),i*cosd(mara),0);
    vt := 0.25*lens*(sind(mara),cosd(mara),0);
    c[coun]1 = vt+vr;
    c[coun]2 = vt+vr+lens*blue;
    c[coun]3 = vr-vt+lens*blue;
    c[coun]4 = vr-vt;
    c[coun]5 = 0.25*(c[coun]1+c[coun]2+c[coun]3+c[coun]4);
    j := i;
  endfor;
  anga := (lens-flen+j)/sray/pif;
  for i=anga-mara step lens/sray/pif until mara:
    coun := incr( coun );
    vr := armb*red + (sray*cosd(i),sray*sind(i),3*sqsd);
    vt := 0.25*lens*(-sind(i),cosd(i),0);
    c[coun]1 = vt+vr;
    c[coun]2 = vt+vr+lens*blue;
    c[coun]3 = vr-vt+lens*blue;
    c[coun]4 = vr-vt;
    c[coun]5 = 0.25*(c[coun]1+c[coun]2+c[coun]3+c[coun]4);
    j := i;
  endfor;
  auxd := lens-pif*(mara-j)*sray;
  vi := armb*red + (sray*cosd(mara),sray*sind(mara),3*sqsd);
  for i=auxd step lens until flen:
    coun := incr( coun );
    vr := vi + (-i*sind(mara),i*cosd(mara),0);
    vt := 0.25*lens*(-sind(mara),cosd(mara),0);
    c[coun]1 = vt+vr;
    c[coun]2 = vt+vr+lens*blue;
    c[coun]3 = vr-vt+lens*blue;
    c[coun]4 = vr-vt;
    c[coun]5 = 0.25*(c[coun]1+c[coun]2+c[coun]3+c[coun]4);
  endfor;
  for i=1 upto coun:
    for j=1 upto 4:
      c[i][j] := eulerrotation( acta, 0, 0, c[i][j] );
    endfor;
    elo[i] = for k=1 upto 4: rp(c[i][k])-- endfor cycle;
    inde := decimal( i );
    comm := "unfill elo[" & inde & "]; draw elo[" & inde & "];";
    getready( comm, c[i]5 );
  endfor;
  doitnow;
  produce_vga_border
endfig;
end;







