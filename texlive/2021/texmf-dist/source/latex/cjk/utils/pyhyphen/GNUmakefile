# Makefile to create hyphenation patterns for pinyin with tone marks and
# pinyin without tone marks.

BASH = bash
CC = cc
DATE = date
GIT = git
ICONV = iconv -f latin1 -t utf8
RM = rm -f
SED = sed # should be GNU sed

GIT_VERSION = $(shell $(GIT) log --format=%H -1 HEAD --)
LEFTHYPHENMIN = $(strip $(shell $(SED) 's/^\(..\).*/\1/;q' < pinyin.tr))
RIGHTHYPHENMIN = $(strip $(shell $(SED) 's/^..\(..\).*/\1/;q' < pinyin.tr))
TODAY = $(shell $(DATE) '+%Y-%m-%d')


all: hyph-zh-latn-pinyin.tex \
     hyph-zh-latn-tonepinyin.tex

pinyin: pinyin.c
	$(CC) -o $@ $<

pinyin.dic: pinyin
	./$< > $@

pattern.rules pinyin.pat: make-patterns.sh \
                          pinyin.dic \
                          pinyin.tr
	$(BASH) make-patterns.sh pinyin.dic pinyin.tr \
        && $(ICONV) < pattern.2 > pinyin.pat

tonepinyin.pat: pinyin.pat
	LANG=en_US.UTF-8 \
        $(SED) -e 's/\(.*\)a\(.*\)/\1a\2 \1ā\2 \1á\2 \1ǎ\2 \1à\2/' \
               -e 's/\(.*\)e\(.*\)/\1e\2 \1ē\2 \1é\2 \1ě\2 \1è\2/' \
               -e 's/\(.*\)i\(.*\)/\1i\2 \1ī\2 \1í\2 \1ǐ\2 \1ì\2/' \
               -e 's/\(.*\)o\(.*\)/\1o\2 \1ō\2 \1ó\2 \1ǒ\2 \1ò\2/' \
               -e 's/\(.*\)u\(.*\)/\1u\2 \1ū\2 \1ú\2 \1ǔ\2 \1ù\2/' \
               -e 's/\(.*\)ü\(.*\)/\1ü\2 \1ǖ\2 \1ǘ\2 \1ǚ\2 \1ǜ\2/' \
               < $< > $@

hyph-zh-latn-pinyin.tex: hyph-zh-latn-pinyin.in \
                         pattern.rules \
                         pinyin.pat
	$(SED) -e 's/@WITH@/without/' \
               -e 's/@DATE@/$(TODAY)/' \
               -e 's/@GIT_VERSION@/$(GIT_VERSION)/' \
               -e 's/@LEFTHYPHENMIN@/$(LEFTHYPHENMIN)/' \
               -e 's/@RIGHTHYPHENMIN@/$(RIGHTHYPHENMIN)/' \
               -e '/used patgen parameters/ r pattern.rules' \
               -e '/\\patterns/ r pinyin.pat' \
               < $< > $@

hyph-zh-latn-tonepinyin.tex: hyph-zh-latn-pinyin.in \
                             pattern.rules \
                             tonepinyin.pat
	$(SED) -e 's/@WITH@/with/' \
               -e 's/@DATE@/$(TODAY)/' \
               -e 's/@GIT_VERSION@/$(GIT_VERSION)/' \
               -e 's/@LEFTHYPHENMIN@/$(LEFTHYPHENMIN)/' \
               -e 's/@RIGHTHYPHENMIN@/$(RIGHTHYPHENMIN)/' \
               -e '/used patgen parameters/ r pattern.rules' \
               -e '/\\patterns/ r tonepinyin.pat' \
               < $< > $@

clean:
	$(RM) pinyin pinyin.dic pinyin.pat tonepinyin.pat
	$(RM) pattern.*
	$(RM) pattmp.*

# Local Variables:
# mode: makefile
# coding: utf-8
# End:
