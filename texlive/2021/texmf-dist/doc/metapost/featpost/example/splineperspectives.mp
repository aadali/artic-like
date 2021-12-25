% splineperspectives.mp
% L. Nobre G.
% Troy Henderson
% Przemek Koprowski
% Manjusha Joshi
% 2012

prologues := 1;
color f, vecU, vecV, Perp;
numeric EPS, iterations;
EPS:=1/25;
iterations = 0;

    vardef initializeperspective( expr Focus ) =
      f := Focus;
      vecV := N( (-Y(f), X(f), 0) );
      vecU := ncrossprod( f, vecV );
      Perp := N( f )
    enddef;

    def X(expr A) = redpart A enddef;
  
    def Y(expr A) = greenpart A enddef;
  
    def Z(expr A) = bluepart A enddef;

    def conorm(expr A) = ( X(A) ++ Y(A) ++ Z(A) ) enddef;

    def N(expr A) =
        begingroup
            save M, exitcolor;
            numeric M;
            color exitcolor;
            M = conorm( A );
            if M > 0:
                exitcolor = ( X(A)/M, Y(A)/M, Z(A)/M );
            else:
                exitcolor := black;
            fi;
            ( exitcolor )
        endgroup
    enddef;

    def cdotprod(expr A, B) = 
        ( X(A)*X(B) + Y(A)*Y(B) + Z(A)*Z(B) )
    enddef;

    def ccrossprod(expr A, B) = 
        ( Y(A)*Z(B) - Z(A)*Y(B), 
          Z(A)*X(B) - X(A)*Z(B), 
          X(A)*Y(B) - Y(A)*X(B) )
    enddef;

% The dotproduct of two normalized vectors is the cosine of the angle 
% they form.

    def ndotprod(expr A, B) = 
        begingroup
            save a, b;
            color a, b;
            a = N(A);
            b = N(B);
            ( ( X(a)*X(b) + Y(a)*Y(b) + Z(a)*Z(b) ) )
        endgroup
    enddef;

% The normalized crossproduct of two vectors. 
% Also check getangle below.

    def ncrossprod(expr A, B) = 
        N( ccrossprod( A, B ) )
    enddef;

    def rp(expr R) =
        begingroup
	  numeric verticalcoordinate, horizontalcoordina;
	  f-vecV*horizontalcoordina-vecU*verticalcoordinate=whatever*(f-R);
	  ( 150*(horizontalcoordina,verticalcoordinate) )
        endgroup
    enddef;

    def line( expr Ang ) =
      begingroup
	save a, b, c;
	numeric a, b, c;
	a = (2-(1 ++ cosd(Ang))*cosd(3*Ang))*cosd(Ang);
	b = (2-(1 ++ cosd(Ang))*cosd(3*Ang))*sind(Ang);
	c =1.5+(1 ++ cosd(Ang))*sind(3*Ang);
	( (a,b,c) )
      endgroup
    enddef;
    
    vardef rationalnobreg(expr A,B,C,D) =
      save P,Q,E,wmin,wMax,we;
      color P[],Q[],E;
      numeric we[];
      we[0]:=cdotprod( Perp, f-A );
      we[1]:=cdotprod( Perp, f-B );
      we[2]:=cdotprod( Perp, f-C );
      we[3]:=cdotprod( Perp, f-D );
      wmin:=min(we[0],we[1],we[2],we[3]);
      wMax:=max(we[0],we[1],we[2],we[3]);
      if ((1-wmin/wMax) > EPS):
	P[0]:=A;
	P[1]:=1/2[A,B];
	E:=1/2[B,C];
	Q[2]:=1/2[C,D];
	Q[3]:=D;
	P[2]:=1/2[P[1],E];
	Q[1]:=1/2[E,Q[2]];
	P[3]:=1/2[P[2],Q[1]];
	Q[0]:=P[3];
	iterations := incr( iterations );
	rationalnobreg(P[0],P[1],P[2],P[3]) & rationalnobreg(Q[0],Q[1],Q[2],Q[3])
      else:
	rp(A) .. controls rp(B) and rp(C) .. rp(D)
      fi
    enddef;
  
    def casteljau( expr Za, Zb, Zc, Zd, Pt ) = %%%%%%%%%%%%%%%%%%% 2D or 3D
      begingroup
	save A, B, C, D;
	numeric A, B, C, D;
	A = (1-Pt)**3;
	B = 3*((1-Pt)**2)*Pt;
	C = 3*(1-Pt)*(Pt**2);
	D = Pt**3;
	( (A*Za+B*Zb+C*Zc+D*Zd) )
      endgroup
    enddef;
    
    def twothr( expr Z ) = ( xpart Z, ypart Z, 0 ) enddef;
    
    def twotwo( expr Z ) = rp( twothr( Z ) ) enddef;

    def xoy( expr Z ) = rp( ( X(Z), Y(Z), 0 ) ) enddef;

    def yoz( expr W ) = rp( ( 0, Y(W), Z(W) ) ) enddef;

    def xoz( expr W ) = rp( ( X(W), 0, Z(W) ) ) enddef;

    def nextthirty( expr Za, Zb, Zc, Zd, Pt ) = %%% input 3D and return 2D
      begingroup
	save A, B, C, D, Tot, P;
	numeric A, B, C, D, Tot;
	A = ((1-Pt)**3)*cdotprod( Perp, f-Za );
	B = 3*((1-Pt)**2)*Pt*cdotprod( Perp, f-Zb );
	C = 3*(1-Pt)*(Pt**2)*cdotprod( Perp, f-Zc );
	D = (Pt**3)*cdotprod( Perp, f-Zd );
	Tot = A+B+C+D;
	( (A*rp(Za)+B*rp(Zb)+C*rp(Zc)+D*rp(Zd))/Tot )
      endgroup
    enddef;

beginfig(1);
%  initializeperspective(0.35*(3,5,2));
  initializeperspective((1.3,1.31,1.32));
  color w[];
  numeric num, i;
  pen pencontrol, penalytic;
  color colcontrol, colorytic;
  pencontrol = pencircle scaled 4pt;
  penalytic = pencircle scaled 2pt;
  colcontrol = black;
  colorytic = blue+green;
  num = 50;
  w1 = (1,0,0);
  w2 = (0,0,1);
  w3 = (0,1,0);
  w4 = (1,1,1);
  w5 = (1,1,0);
  w6 = (1,0,1);
  w7 = (0,1,1);
  for i=1 upto 3: draw rp(black)--rp(w[i]); endfor;
  draw rp(w1)--rp(w2)--rp(w3)--rp(w4)--rp(w5)--rp(w1)--rp(w6)--
       rp(w2)--rp(w7)--rp(w3)--rp(w5) dashed withdots;
  draw rp(w6)--rp(w4)--rp(w7) dashed withdots;
  draw xoy(w1) for i=1 upto num: ..xoy(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw xoz(w1) for i=1 upto num: ..xoz(casteljau(w1,w2,w3,w4,i/num)) endfor;
  draw rp(w1) for i=1 upto num: ..nextthirty(w1,w2,w3,w4,i/num) endfor
                                  withpen pencontrol withcolor colcontrol;
  draw rp(w1) for i=1 upto num: ..rp(casteljau(w1,w2,w3,w4,i/num)) endfor
                                  withpen penalytic withcolor colorytic;
  draw rationalnobreg(w1,w2,w3,w4);
  show iterations;
  draw nextthirty(w1,w2,w3,w4,0.5) withpen pencontrol withcolor red;	      
endfig;

beginfig(2);
  initializeperspective((3,5,4));
  color node[], pre[], pos[], a, b, c;
  numeric i, j, k, l;
  path td[];
  pen penmark;
  pair aux[];
  penmark = pencircle scaled 2mm;
  j = 0;
  draw for i=1 step 3 until 360: xoy(line(i)).. endfor cycle;
  draw for i=1 step 3 until 360: rp(line(i)).. endfor cycle;
  for i=10 step 23 until 360:
    j := incr( j );
    node[j] = line( i );
    draw rp(node[j]) withpen penmark;
  endfor;
  td1 = for i=1 upto j: (X(node[i]),Y(node[i])).. endfor cycle;
  td2 = for i=1 upto j: (X(node[i]),Z(node[i])).. endfor cycle;
  td3 = for i=1 upto j: (Y(node[i]),Z(node[i])).. endfor cycle;
  for i=1 upto j:
    l := i-1;
    for k=1 upto 3:
      aux[k] := precontrol l of td[k];
    endfor;
    a := (xpart (aux1 - point l of td1),ypart (aux1 - point l of td1),0);
    b := (xpart (aux2 - point l of td2),0,ypart (aux2 - point l of td2));
    c := (0,xpart (aux3 - point l of td3),ypart (aux3 - point l of td3));
    pre[i] = node[i]+0.5*(a+b+c);
    draw rp(pre[i]) withpen penmark withcolor blue;
    for k=1 upto 3:
      aux[k] := postcontrol l of td[k];
    endfor;
    a := (xpart (aux1 - point l of td1),ypart (aux1 - point l of td1),0);
    b := (xpart (aux2 - point l of td2),0,ypart (aux2 - point l of td2));
    c := (0,xpart (aux3 - point l of td3),ypart (aux3 - point l of td3));
    pos[i] = node[i]+0.5*(a+b+c);    
    draw rp(pos[i]) withpen penmark withcolor red;
    draw rp(pre[i])--rp(pos[i]) withcolor green;
  endfor;
endfig;

beginfig(3);
  initializeperspective((7,11,9));
  color node[], pre[], pos[];
  numeric i, j, param, k, num;
  pen penmark;
  param = 0.333;
  num = 12;
  penmark = pencircle scaled 2mm;
  j = 1;
  draw for i=1 step 3 until 360: xoy(line(i)).. endfor cycle;
  draw for i=1 step 3 until 360: rp(line(i)).. endfor cycle;
  for i=10 step 23 until 360:
    j := incr( j );
    node[j] = line( i );
    draw rp(node[j]) withpen penmark;
  endfor;
  node[1] = node[j];
  node[0] = node[j-1];
  node[j+1] = node[2];
  for i=2 upto j:
    %
    pre[i] = node[i]-param*conorm(node[i]-node[i-1])*N(node[i+1]-node[i-1]);
    %
  endfor;
  pre[1] = pre[j];
  for i=1 upto j-1:
    %
    pos[i] = node[i]+param*conorm(node[i+1]-node[i])*N(node[i+1]-node[i-1]);
    %
  endfor;
  pos[j] = pos[1];
  for i=1 upto j-1: 
    draw rp(pre[i]) withpen penmark withcolor blue;
    draw rp(pos[i]) withpen penmark withcolor red;
    draw rp(pre[i])--rp(pos[i]) withcolor green;
    draw
      for k=0 upto num-1:
        rp(casteljau(node[i],pos[i],pre[i+1],node[i+1],k/num))..
      endfor rp(node[i+1]) withcolor (red+blue);
  endfor;
endfig;

end.
