% signalvertexSD.mp
% L. Nobre G. 
% 2003

input featpost3Dplus2D;

% Draw a torus with NP dots, a section of 3.14*RS^2
% and a total diameter of 2*(RB+RS). The dots are placed 
% on a single closed line that turns around the $z$ axis
% NB times and around the section NS times. 
  
f := 0.45f;
                
%SphericalDistortion := true;  Spread := 11;

beginfig(1);
    NP:=2000;
    RB:=1;
    RS:=0.2;
    R3:=0.0;
    NB:=8;
    NS:=5;
    jB:=360*NB/NP;
    jS:=360*NS/NP;
    for i=1 upto NP:
        phi   :=i*jS;
        theta :=i*jB;
        first :=(RB+(RS+R3*sind(4*theta))*cosd(phi))*cosd(theta);
        second:=(RB+(RS+R3*sind(4*theta))*cosd(phi))*sind(theta);
        third :=    (RS+R3*sind(4*theta))*sind(phi);
        signalvertex( (first,second,third), 1, black );
    endfor;
endfig;

end;





