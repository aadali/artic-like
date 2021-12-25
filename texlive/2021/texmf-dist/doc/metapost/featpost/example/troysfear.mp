% troysfear.mp 
% L. Nobre G.
% 2008

input featpost3Dplus2D;

% verbatimtex 
% \documentclass{article} 
% \usepackage{beton} 
% \usepackage{concmath} 
% \usepackage{ccfonts} 
% \begin{document} 
% etex

prologues := 1;

beginfig(1);
  f := (10,-5,5.44);
  Spread := 120;
  viewcentr := black;
  numeric i, numl, aux, thickf;
  color poincol, gammacnt, raycol;
  numl = 20;
  thickf = 0.05;
  poincol = blue;
  raycol = (1,1,0.1);
  color poin[];
  poin1 = ( 0.0000000000000000, 0.0000000000000000, 1.0000000000000000);
  poin2 = ( 0.7333266247606663, 0.0000000000000000, 0.6798765045338225);
  poin3 = ( -0.4170317493123703, -0.8933183244372310, -0.1675317620337181);
  poin4 = ( -0.2314815726239148, 0.6208617150619569, 0.7489639592836793);
  poin5 = ( 0.1537652646366442, -0.9604012907848274, 0.2323480235553130);
  poin6 = ( -0.3863954912627506, -0.6607125590945075, 0.6435506495891403);
  poin7 = ( -0.7710058807464431, 0.5605065386692130, 0.3022951404893231);
  poin8 = ( -0.0959044889866723, 0.9868269446283826, 0.1302877981532283);
  poin9 = ( 0.0997534549176198, 0.7727171212494309, -0.6268632217318151);
  poin10 = ( 0.8172909883532864, -0.5598423938515633, 0.1364255635980997);
  poin11 = ( 0.2373103227489228, 0.0812688716483729, -0.9680285022755282);
  poin12 = ( 0.3405602349795320, -0.8413398111607557, -0.4197213939110916);
  poin13 = ( -0.9610728898061991, 0.1533664531864700, -0.2298208683226212);
  poin14 = ( -0.8870795087724005, -0.4278302969039353, 0.1733527679825362);
  poin15 = ( 0.7212726731262007, 0.4012160791697596, -0.5646161428944713);
  poin16 = ( -0.5528113739219546, 0.7452763339754855, -0.3727771061620893);
  poin17 = ( 0.5797275131891773, 0.8147563131365292, -0.0093893905905605);
  poin18 = ( -0.0762908334114027, -0.5230879059144918, -0.8488573210047642);
  poin19 = ( 0.4559835073369932, 0.6170900853351784, 0.6413102740621552);
  poin20 = ( 0.7633512867021243, -0.3106468355352665, -0.5663862256995101);
  poin21 = ( -0.4508345799824528, 0.2265357804581852, -0.8633827202719813);
  poin22 = ( -0.7074496211764032, -0.3987032948435680, -0.5835672336400034);
  poin23 = ( 0.9849663581348018, 0.1640719469282586, 0.0540524705618168);
  poin24 = ( 0.3404137309548133, -0.5605840215601572, 0.7548934007850813);
  poin25 = ( -0.6907194772069863, -0.0072100091830313, 0.7230868686226352);
  for i=1 upto 25:
%    draw rp( N( poin[i] ) )--rp( poin[i] ) withcolor raycol;
%    signalvertex( poin[i], thickf, poincol );
    fill rigorouscircle( poin[i], poin[i], thickf ) withcolor poincol;
  endfor;
  tropicalglobe( numl, black, 1, blue );
  for i=1 upto numl:				% longitudes
    aux		:= (i-1)*180/numl;
    gammacnt 	:= (cosd(aux),sind(aux),0);
    draw spatialhalfcircle(black,gammacnt,1,true);
  endfor;
  for i=1 upto 25:
    if ( conorm( f-poin[i] ) < conorm( f ) ) and ( conorm( poin[i] ) > 0.99 ):
%      draw rp( N( poin[i] ) )--rp( poin[i] ) withcolor raycol;
%      signalvertex( poin[i], thickf, poincol );
    fill rigorouscircle( poin[i], poin[i], thickf ) withcolor poincol;
    fi;
  endfor;
endfig;

% verbatimtex 
% \end{document} 
% etex

end;
