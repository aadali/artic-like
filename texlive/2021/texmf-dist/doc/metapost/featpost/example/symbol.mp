% symbol.mp
% G. Morais
% 2008

verbatimtex 
\documentclass{article}
\usepackage{palatino}
\begin{document} 
etex

prologues := 3;

beginfig(1)

  numeric u ; u = 110 ;
  numeric int ; int = .053u ;
  numeric ext ; ext = .1u ;

  color mycolor ;
  mycolor := (.45 , .45 , .45 ) ;

  color intcontorno ;
  intcontorno := (.23 ,.23 ,.23 ) ;
  
  color ground ;
  ground := ( .09 , .19 , .25 ) ;

  color sky ;
  sky := ( .24 , .30 , .86 ) ;
  
  path contorno[] ; 
  contorno[1] := (-2,2)*u..(0,-1)*u..(2,3)*u ;

  z0 = point 0 of contorno[1] ;
  z1 = point 1.8 of contorno[1] ;
  z2 = point 0.5 of contorno[1] ;
  z3 = point 2 of contorno[1] ;

  contorno[2] := z0{(1,2)}..z3 ;

  path interior[] ;
  interior[1] := z1..(1,0)*u..z2 ;
  z4 = point 1.5 of interior[1] ;

  interior[2] := z4{right}..{right}z3 ;
  z5 = point 0.4 of interior[2] ;

  interior[3] := z0{(1,2)}..(.5,.35)*u..{(1,2)}z5 ;

  path auxiliar[] ;
  auxiliar[1] := (2,3)*u .. (0,-1)*u .. (-2,2)*u ;

  auxiliar[2] := interior[2] .. subpath(2,1.8) of contorno[1] ..
                 interior[1] cutafter z4 .. cycle ;

  path ciclo[] ;

  ciclo[1] := interior[1] .. ( contorno[1] cutbefore z2 ) cutafter z1 ..cycle ;
  ciclo[2] := interior[2] .. subpath(2,1.8) of contorno[1] ..
              interior[1] cutafter z4 .. cycle ;
  ciclo[3] := buildcycle ( interior[3] , interior[2] , interior[1] , contorno[1] ) ;
  ciclo[4] := buildcycle ( interior[3] , interior[2] , contorno[2] ) ;

  fill ciclo[1] withcolor ground ;
  fill ciclo[2] withcolor .2blue ;
  fill ciclo[3] withcolor .3red ;
  fill ciclo[4] withcolor sky ;
  
  pickup pencircle scaled int ;
  draw interior[1] withcolor intcontorno ;
  draw interior[2] withcolor intcontorno ;
  draw interior[3] withcolor intcontorno ;

  pickup pencircle scaled ext ;
  draw contorno[1] withcolor mycolor ;
  draw contorno[2] withcolor mycolor ;

  color letras[] ;
  letras[1] := ( .59 , .11 , .06 ) ;
  letras[2] := ( .07 , .02 , .43 ) ;
  draw btex GNUPOST etex
  xscaled .112u yscaled .06u
  shifted ( -2.8u , -1.6u )
  withcolor letras[1] ;
  draw btex WE SHARE OUR VISION etex scaled .052u
  shifted ( -2.75u , -2.05u )
  withcolor letras[2] ;

%  setbounds currentpicture to boundingbox currentpicture enlarged .05u ;
  
endfig;

verbatimtex 
\end{document} 
etex

end ;