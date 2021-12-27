% keplernewtonfeynman.mp
% L. Nobre G. 
% 2014

%input featpost3Dplus2D;

prologues := 0;

%verbatimtex \documentclass{article}\usepackage{mathpazo}\begin{document} etex

ahangle := 30;

beginfig(1);
  numeric u, abratio, foc, pang, auxang, refdist, stepang, factor, thang;
  pair secofod, ocirpoi, midpoi, tpoi, dirplu, dirmin;
  pen poipen;
  color mygrey;
  path thecircle;
  thang = -45;
  stepang = 4;
  poipen = pencircle scaled 1mm;
  mygrey = 0.5*(red+green);
  u = 1cm;
  %abratio = 0.5*(1+sqrt(5));
  abratio = 1.2;
  foc = u*(abratio +-+ 1);
  refdist = 2*u*abratio;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  factor = 1.02;
  secofod := 2*foc*dir(thang);
  thecircle = fullcircle scaled 2refdist;
  pickup pencircle scaled 0;
  %draw fullcircle scaled 2u xscaled abratio shifted (foc,0);
  %draw thecircle;
  draw secofod withpen poipen withcolor blue;
  for pang = stepang step stepang until 360:
    ocirpoi := refdist*dir(pang);
    midpoi := 0.5[secofod,ocirpoi];
    auxang := angle(ocirpoi-secofod);
    dirplu := dir(auxang+90);
    dirmin := dir(auxang-90);
    pair tpoi;
    tpoi = whatever*dir(pang);
    tpoi = whatever[midpoi,midpoi+dirplu];
    %draw ocirpoi--secofod withcolor mygrey;
    %draw tpoi withpen poipen withcolor red;
    drawarrow tpoi--(tpoi+abs(midpoi-secofod)*dirplu/factor);
  endfor;
  clip currentpicture to (thecircle scaled (factor**2));
  draw origin withpen poipen withcolor green;
endfig;

verbatimtex \end{document} etex

end.

Newton showed this construction in Book 1, Section 4, Lemma 15, of
Principia. On March 13, 1964, Feynman resurrected the construction and
used it in a lecture, "The Motion of Planets Around the Sun". The
lecture is detailed in a book with audio CD, Feynman's Lost Lecture,
by David and Judith Goodstein. In the lecture, Feynman used the
diagram and differential geometry to prove the planetary laws of
motion.  
In the construction, the green dot is the primary focus of the ellipse
about which the planet orbits; the blue dot is the secondary
focus. The black dot is on a circle at a distance in radius equal in
length to the major axis of an ellipse. A line is drawn from the blue
dot to the black dot and its perpendicular bisector is
constructed. The point where this perpendicular bisector intersects
the line from the green dot (primary focus) to the black dot (circle)
is a point on the ellipse. The perpendicular bisector is tangent to
the ellipse at this point. A further interesting point of the
construction is that the length of the line from the blue dot
(secondary focus) to the black dot (circle) is proportional to the
velocity of the orbiting planet at this point. In the Demonstration,
half this length is represented by the black vector traveling in the
direction of the planet.  
All of the possible ellipses with the given major axis are contained
in the circle. You can adjust the eccentricity and rotation of the
ellipse.  