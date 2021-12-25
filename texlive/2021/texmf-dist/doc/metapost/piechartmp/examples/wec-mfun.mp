%H
%H piechartMP example - world energy consumption 2000
%H (extended MetaFun version)
%H Jens-Uwe Morawski, 2002
%H


%I this MetaPost example requires the MetaFun package
%I (see piechartMP documentation)
%I The purpose of this example is to show how pie-charts
%I can be placed in a figure using the pc_Centre variable.
%I wec-mfun.pdf in the same directory contains the final
%I figure.


%D This example uses the 'externalfigure' macro from
%D the MetaFun package. If piechartMP is loaded after
%D MetaFun then piechertMP detects MetaFun and will
%D draw segment outsides as gradients.
%D Since here many charts will be drawn it can cause
%D an overflow of the MetaPost buffer-size,
%D therefore piechertMP is loaded prior MetaFun.
%D If the buffer-size of your MetaPost installation is
%D large enough then you can load piechartMP after
%D MetaFun and segment outsides will be drawn as gradients.
input piechartmp ;
input metafun ;

%D Helvetica-Bold
% This should work on most systems. If not, comment the next line.
defaultfont := "phvb8r" ;
defaultscale := 0.7 ;

%D 
Scale := 5mm ;

%D define some colors
color oil_color, gas_color, coal_color, nuclear_color, water_color ;
oil_color     := (0.4, 0  , 0  ) ;
gas_color     := (0  , 0.5, 0  ) ;
coal_color    := (0.5, 0.3, 0  ) ;
nuclear_color := (0.8, 0.8, 0.2) ;
water_color   := (0.3, 0.3, 1.0) ;


%D define some points to place the charts
pair place[] ;
place[1] := (2cm,5cm) ;% North America
place[2] := (7.5cm,3.5cm) ;% Asia
place[3] := (5cm,5.3cm) ;% Europe
place[4] := (6.5cm,5.8cm) ;% CIS
place[5] := (5.6cm,4.2cm) ;% Middle East
place[6] := (2.5cm,2.7cm) ;% South America
place[7] := (5cm,3cm) ;% Africa


%D
def AChart (expr Size, Name, Place, Water, Nuclear, Coal, Gas, Oil) =
	begingroup;
	save label_pnt, Radius ;
	pair label_pnt ; numeric Radius ;

	% the 'Size' is the absolute energy consumption
	% that means the 'Size' is proportional to the area
	% of the pie-chart and not the radius.
	Radius := (sqrt Size)*Scale ;

	ResetSegments ; % documented in piechartMP Internals 

		if Water>0:
			Segment(Water, "", water_color) ;
		fi;
		if Nuclear>0:
			Segment(Nuclear, "", nuclear_color) ;
		fi;
		if Coal>0:
			Segment(Coal, "", coal_color) ;
		fi;
		if Gas>0:
			Segment(Gas, "", gas_color) ;
		fi;
		if Oil>0:
			Segment(Oil, "", oil_color) ;
		fi;

		pc_Centre := place[Place] ; % documented in piechartMP Internals 
		PieChart(Radius, 0.3, 65, 90, 0) ;
		label_pnt := SegmentPoint(1,(1,0));
		label.top (Name & "; " & (decimal Size), label_pnt) withcolor white;

	endgroup;
enddef;


beginfig(1);

%D 'externalfigure' is a macro from the MetaFun package.
%D It embeds the external JPG-Background, the file worldmap.jpg
%D should be available from the same location as this
%D MP file.
externalfigure "worldmap.jpg" xscaled 10cm  yscaled 7cm ;

AChart (3.8, "North-America", 1, 2, 6, 25, 32, 35 ) ;
AChart (3.4, "Asia/Oceania", 2, 2, 5, 46, 10, 38 ) ;
AChart (2.6, "Europe", 3, 4, 16, 18, 20, 42 ) ;
AChart (1.3, "CIS", 4, 2, 6, 18, 50, 24 ) ;
AChart (0.5, "Middle East", 5, 0, 0, 0, 35, 65 ) ;
AChart (0.5, "South-America",6 , 15, 0, 5, 15, 65 ) ;
AChart (0.4, "Africa",7 , 0, 0, 37, 18, 45 ) ;


%D the legend
x := 3mm ;
y := 2.5cm ;
fill fullcircle scaled 2mm shifted (x,y) withcolor oil_color;
label.rt ("Oil", (x+2mm,y)) ;
fill fullcircle scaled 2mm shifted (x,(y-5mm)) withcolor gas_color;
label.rt ("Gas", (x+2mm,y-5mm)) ;
fill fullcircle scaled 2mm shifted (x,(y-1cm)) withcolor coal_color;
label.rt ("Coal", (x+2mm,y-1cm)) ;
fill fullcircle scaled 2mm shifted (x,(y-1.5cm)) withcolor nuclear_color;
label.rt ("Nuclear Power", (x+2mm,y-1.5cm)) ;
fill fullcircle scaled 2mm shifted (x,(y-2cm)) withcolor water_color;
label.rt ("Water Power", (x+2mm,y-2cm)) ;


%D some labels
label.top ("(in billion t MCU)", (7cm,0)) ;
defaultfont:="ptmbi8r" ; defaultscale:= 1.2;
label("World Energy Consumption 2000", origin)
	rotated 90 shifted (9.5cm,3cm) withcolor white;


%D Next follows the grid I have used to determine the co-ordinates
%D saved in place[]. Now it is not drawn because of the if..fi
%D statement.
if false:
for i=1 upto 10:
	draw (origin--(0,7cm)) shifted (i*1cm,0) withcolor white ;
endfor;
for i=1 upto 7:
	draw (origin--(10cm,0)) shifted (0,i*1cm) withcolor white ;
endfor;
fi;
endfig;
end