%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ten plik jest czê¶ci± pakietu MPATTERN
%
% Autor: Piotr Bolek
% wersja: 0.5 (Jun 25, 2001)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Opis
====

Pakiet MPATTERN s³u¿y do definiowania i u¿ywania wype³nieñ wzorkowych
w programie MetaPost. W pakiecie wykorzystywana jest dostêpna w jêzyku
PostScript Level 2 Wzorkowa Przestrzeñ Kolorów (Pattern Color Space).

Zawarto¶æ:
- mpattern.mp: makra dla programu MetaPost s³u¿±ce do  definiowania
  i u¿ywania wzorków.
- README: krótki opis po angielsku.
- README.PL: krótki opis po polsku.
- test.mp, test.tex, test.ps: przyk³ady u¿ycia.

Sposób u¿ycia
=============

W pliku ¼ród³owym nale¿y na pocz±tku w³±czyæ poleceniem input plik
mpattern.mp, a nastêpnie przetworzyæ plik standardowym programem
mpost. 

Interfejs u¿ytkownika
=====================

Definiowanie wzorku
--------------------

Wzorki s± definiowane przy u¿yciu konstrukcji:

beginpattern(<nazwa_wzorku>);
...
endpattern;

Parametr <nazwa_wzorku> jest nazw±, która s³u¿y do identyfikacji
wzorku. Nazwa ta mo¿e zawieraæ litery i znaki podkre¶lenia ,,_'', ale
nie mo¿e zawieraæ cyfr. NIE nale¿y podawaæ tej nazwy w cudzys³owach.

W definicji wzorku (w dowolnym miejscu miêdzy beginpattern a
endpattern) mo¿na podaæ dodatkowe polecenia: 
- patternbbox(<pair>,<pair>)
- patternbbox((<numeric>,<numeric>,<numeric>,<numeric>) 
  Jawne podanie wspó³rzêdnych prostok±ta ograniczaj±cego komórkê wzorka

- patternxstep(<numeric>)
- patternystep(<numeric>)
- patternstep(<numeric>,<numeric>)
- patternstep(<pair>)
  Odstêpy miêdzy komórkami wzorka

- patterntransform(<transform>)
  Modyfikacja kszta³tu komórki wzorka

- patterncolor(0<= <numeric> <= 1)
  Jawne podanie koloru zwi±zanego z danym wzorkiem

Ka¿dy wzorek jest tymczasowo zapisywany jako rysunek o du¿ym numerze
(domy¶lnie 999), ale je¶li u¿ytkownik wykorzysta ten numer do swojego
rysunku, to dla wzorku bêdzie u¿yty najwiêkszy wolny numer mniejszy od
999.

U¿ywanie wzorku
---------------

Do wykorzystania zdefiniowanego wzorku s³u¿y makro dzia³aj±ce jak
operator dwuargumentowy (primary) withpattern. Sposób u¿ycie jest
nastêpuj±cy:

<¶cie¿ka> withpattern <nazwa_wzorku>;

gdzie <¶cie¿ka> jest dowolnym wyra¿eniem okre¶laj±cym zamkniêt±
¶cie¿kê, a <nazwa_wzorku> jest nazw± podan± w definicji wzorku.

Przyk³ad
--------

Oto kod s³u¿±cy do narysowania dwóch kó³ek wype³nionych ró¿nymi
szachownicami.

input mpattern;

beginpattern(checker);
 d:=4mm;
 fill unitsquare scaled d;
 fill unitsquare scaled d shifted (d,d);
 patternbbox(0,0,2d,2d);
endpattern;

beginpattern(rotated_checker);
 fill unitsquare scaled 4mm rotated 45 withcolor .7white;
endpattern;

beginfig(11);
 path p;
 z1=(1cm,1cm);
 p=fullcircle scaled 3cm;
 fill p withpattern checker;
 unfill p shifted z1;
 fill p shifted z1 withpattern rotated_checker;
endfig;

Zasady rozpowszechniania
========================

Pakiet MPATTERN jest dobrem wspólnym (ang. ,,public domain''). W
przypadku modyfikacji lub redystrybucji pakietu mo¿na stosowaæ dowoln±
aktualn± licencjê GPL.

Problemy
========

06-2001: Mog± byæ problemy z wype³nianiem wzorkami napisów. Aktualne
	 rozwi±zanie jest prowizoryczne.

Zmiany
======

06-2001: Dziêki inspiracji tutorialem Piotra Strzelczyka na
	 konferencji BachoTeX'2001 wyeliminowana zosta³a konieczno¶æ
	 postprocesingu z wykorzystaniem Perla. Podmiana kolorów na
	 wzorki jest obecnie realizowana z wykorzystaniem s³owników na
	 poziomie PostScriptu.

Kontakt z autorem
=================

e-mail: P.Bolek@gust.org.pl
