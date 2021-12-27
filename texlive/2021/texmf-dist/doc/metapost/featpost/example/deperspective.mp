% deperspective.mp
% L. Nobre G.
% 2006

input featpost3Dplus2D;

defaultscale := 1.5;

verbatimtex \documentclass{article}\begin{document} etex

numeric cont, cmax;
cmax = 4;
color focus[];
focus1 =     (5,  3.5,1);
focus2 = 0.4*(5,  3.5,3);
focus3 = 0.4*(5,  1.5,1.5);
focus4 = 0.4*(2.5,5,  1.5);
for cont= 1 upto cmax:
  f := focus[cont];
  beginfig(cont);
    numeric gridnum, axesmarg, xlen, ylen;
    gridnum = 15;
    axesmarg = 0.5;
    xlen = 1.0;
    ylen = 1.5;
    cartaxes(xlen+axesmarg,ylen+axesmarg,1);
    color va, vb, vc, vd;
    
    va = (   0,   0,0);
    vb = (   0,ylen,0);
    vc = (xlen,ylen,0);
    vd = (xlen,   0,0);
    
    z3 = rp(va);
    z4 = rp(vb);
    z5 = rp(vc);
    z6 = rp(vd);
    
    draw z3--z4--z5--z6--cycle dashed evenly;
    z1 = whatever[z3,z4];
    z1 = whatever[z6,z5];
    z2 = whatever[z6,z3];
    z2 = whatever[z5,z4];
    draw z1--z3 dashed evenly;
    draw z1--z6 dashed evenly;
    draw z2--z3 dashed evenly;
    draw z2--z4 dashed evenly;
    
    numeric aonei, atwoi, bonei, btwoi, conei, ctwoi, donei, dtwoi;
    aonei = 1/abs(z3-z1);
    atwoi = 1/abs(z3-z2);
    bonei = 1/abs(z4-z1);
    btwoi = 1/abs(z4-z2);
    conei = 1/abs(z5-z1);
    ctwoi = 1/abs(z5-z2);
    donei = 1/abs(z6-z1);
    dtwoi = 1/abs(z6-z2);
    
    numeric aone, atwo, bone, btwo, cone, ctwo, done, dtwo;
    aone = abs(z3-z1);
    atwo = abs(z3-z2);
    bone = abs(z4-z1);
    btwo = abs(z4-z2);
    cone = abs(z5-z1);
    ctwo = abs(z5-z2);
    done = abs(z6-z1);
    dtwo = abs(z6-z2);

    numeric i, wei, aux, good;
    for i=0 upto gridnum-1:
      wei := (0.5+i)/gridnum;
      draw z1--rp(wei[vb,vc]) withcolor red;
      draw z2--rp(wei[vd,vc]) withcolor red;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The following lines
%%%%%%%%%%%%%%%%%%%%%% show how rounding error can destroy accuracy!
%    aux := dtwoi/((1-wei)*atwoi+wei*dtwoi);
%    good := dtwoi/(atwoi+wei*(dtwoi-atwoi));
%    show (aux,good);
%    draw z1--(good[z2,z6]);
%    draw z2--((bonei/((1-wei)*aonei+wei*bonei))[z1,z4]);
      draw z1--((1/((1-wei)*abs(z6-z2)/abs(z3-z2)+wei))[z2,z6]) withcolor blue;
      draw z2--((1/((1-wei)*abs(z4-z1)/abs(z3-z1)+wei))[z1,z4]) withcolor blue;
    endfor;
    labels.top(1,2,3,4,5,6);
    show ((atwo/dtwo)*(ctwo/btwo),(aone/bone)*(cone/done));
%show ylen/xlen-(bonei-aonei)/(dtwoi-atwoi);
%show (abs(z2-z6)/abs(z1-z4))*(abs(z1-z4)-abs(z1-z3))/(abs(z2-z6)-abs(z2-z3));
%show sqrt(abs(z3-z4))*sqrt(abs(z5-z6))/(sqrt(abs(z6-z3))*sqrt(abs(z5-z4)));
  endfig;
endfor;
  
verbatimtex \end{document} etex

end;
