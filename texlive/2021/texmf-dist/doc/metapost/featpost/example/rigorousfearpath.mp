% rigorousfearpath.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

beginfig(1);
    numeric sfr, radfac, ind;
    color xv, yv, zv, Ve[];
    path p[];
    f := 2*(5,2,3);
    Spread := 50;
    csr := 0.47;
    clr := 0.57;
    xv := (csr,0,0);
    yv := (0,csr,0);
    zv := (0,0,csr);
    Ve1 := (-1,-1,-1);
    Ve2 := (-1,1,-1);
    Ve3 := (1,-1,-1);
    Ve4 := (1,1,-1);
    Ve5 := (-1,-1,1);
    Ve6 := (-1,1,1);
    Ve7 := (1,-1,1);
    Ve8 := (1,1,1);
    Ve9 := (0,0,0);
    for ind=1 upto 8:
	p[ind] := rigorousfearpath( Ve[ind], csr );
    endfor;
    p9 := rigorousfearpath( Ve9, clr );
    unfill p1;
    draw p1;
    draw rp(Ve1+yv)--rp(Ve2-yv);
    draw rp(Ve1+xv)--rp(Ve3-xv);
    unfill p2;
    draw p2;
    unfill p3;
    draw p3;
    draw rp(Ve4-yv)--rp(Ve3+yv);
    draw rp(Ve4-xv)--rp(Ve2+xv);
    unfill p4;
    draw p4;
    draw rp(Ve1+zv)--rp(Ve5-zv);
    draw rp(Ve2+zv)--rp(Ve6-zv);
    unfill p9;
    draw p9;
    draw rp(Ve3+zv)--rp(Ve7-zv);
    draw rp(Ve4+zv)--rp(Ve8-zv);
    unfill p5;
    draw p5;
    draw rp(Ve5+yv)--rp(Ve6-yv);
    draw rp(Ve5+xv)--rp(Ve7-xv);
    unfill p6;
    draw p6;
    unfill p7;
    draw p7;
    draw rp(Ve8-yv)--rp(Ve7+yv);
    draw rp(Ve8-xv)--rp(Ve6+xv);
    unfill p8;
    draw p8;
endfig;

verbatimtex \end{document} etex

end;

