DVIINCL: made in BOP s.c., bop@bop.com.pl; public domain software
[This is `README' in Polish (slash notation)]

Chcia/loby si/e czasem w/l/aczy/c do obrazka METAPOST-owego
stron/e z przygotowanego wcze/sniej pliku .dvi (zw/laszcza
w kontek/scie prezentacji Bachotkowo-PDF-owo-typograficznej).
Czasem mo/ze to by/c prostsze ni/z wciskanie kodu
TeX-owego mi/edzy znaczniki `btex' i `etex',
szczeg/olnie w przypadku plik/ow wielostronicowych.

Marcin Woli/nski wpad/l na pomys/l, by do tego celu wykorzystywa/c
program DVITOMP. Problem wszak/ze w tym, /ze DVITOMP pstrzy
plik wynikowy znacznikami `mpxbreak', nie tolerowanymi
przez METAPOST-a w normalnym trybie przetwarzania.
Ale nie ma /zadnych przeciwwskaza/n, /zeby wynik programu
DVITOMP wst/epnie przetworzy/c, cho/cby AWK-iem lub PERL-em,
a najlepiej... METAPOST-em!

BOP proponuje nast/epuj/acy spos/ob post/epowania:

1. Plik .dvi konwertujemy programem DVITOMP:
   dvitomp plik_DVI plik_MPX

2. Przetwarzamy plik MPX za pomoc/a METAPOST-a:
   mpost \vardef ifi = str plik_MPX enddef; vardef ofi = str plik_MPY enddef; input fix_brks.mp
   (por. za/l/aczony testit.bat dla DOS/WIN32 lub skrypt testit.sh dla U*X);
   Powstaje plik MPY, chytrze zmieniony w stosunku do pliku MPX -- wstawione
   zostaj/a po prostu instrukcje warunkowe.

   Program `fix_brks.mp' korzysta z za/lo/zenia, /ze operacja `mpxbreak'
   umieszczana jest bez zb/ednych spacji w osobnym wierszu. Gdyby jaki/s
   wariant DVITOMP-a dzia/la/l inaczej, trzeba by zmodyfikowa/c `fix_brks.mp'.

3. W pliku METAPOST-owym, w kt/orym chcemy stron/e (strony) wykorzysta/c,
   nale/zy skorzysta/c z makierka `fix_page' zdefiniowanego
   w pliczku `fix_page.mp', np. tak:
   ...
   picture p;
   p:=fix_page(10,"plik.mpy"); % ,10' to przyk/ladowy numer strony,
                               % "plik.mpy" to przyk/ladowa nazwa
   draw p;
   ...

4. I ju/z.

Komentarze mile widziane.

Kontakt: bop@bop.com.pl
