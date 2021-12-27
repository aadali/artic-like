% tdarrow.mp
% L. Nobre G.
% 2004

input featpost3Dplus2D;

beginfig(1);
  numeric size, num, i, j, coos[];
  string comm;
  color tmpc, fro, tip;;
  size = 1.5;
  num = 250;
  draw rp(black)--rp((size,0,0))--rp((size,size,0))--rp((0,size,0))--cycle;
  for i=1 upto num:
    for j=1 upto 6:
      coos[j] := uniformdeviate(size);
    endfor;
    fro := (coos1,coos2,coos3);
    tip := (coos4,coos5,coos6);
    comm := "tdarrow(" & cstr(fro) & "," & cstr(tip) & ");";
    tmpc := 0.5[fro,tip];
    getready( comm, tmpc );
  endfor;
  doitnow;
endfig;
end.
