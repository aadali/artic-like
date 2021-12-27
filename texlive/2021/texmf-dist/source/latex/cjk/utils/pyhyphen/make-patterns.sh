#!/bin/bash
# -*- coding: utf-8 -*-

#
# This script generates patterns for pinyin.
#
# Usage:
#
#   sh make-patterns.sh pinyin.dic pinyin.tr
#
#
# Input: pinyin.dic  List of all two-syllable combinations of pinyin
#                    syllables (which don't need an apostrophe), with
#                    hyphens inbetween.
#        pinyin.tr   Translation file for patgen.
#
# Output: pattmp.[1-2]       patgen results.
#         pattern.[0-2]      patterns -- pattern.2 is the final pattern file.
#         pattern.[0-2].log  Log files.
#         pattern.rules      A compact description of the used patgen
#                            parameters.
#


# patgen parameters for levels one and two.

hyph_start_finish[1]='1 1'
hyph_start_finish[2]='2 2'

pat_start_finish[1]='2 5'
pat_start_finish[2]='2 5'

good_bad_thres[1]='1 1 1'
good_bad_thres[2]='1 2 1'


# Create empty start patterns.
rm -f pattern.0
touch pattern.0

printf "%%\n" > pattern.rules

for i in 1 2; do

  # Create patterns for current level.  Control instructions for patgen are
  # sent using a pipe.
  printf "%s\n%s\n%s\n%s" "${hyph_start_finish[$i]}" \
                          "${pat_start_finish[$i]}" \
                          "${good_bad_thres[$i]}" \
                          "y" \
  | patgen $1 pattern.$(($i-1)) pattern.$i $2 \
  | tee pattern.$i.log

  # Collect the used patgen parameters.
  printf "%%   %s | %s | %s\n" "${hyph_start_finish[$i]}" \
                               "${pat_start_finish[$i]}" \
                               "${good_bad_thres[$i]}" \
  >> pattern.rules

done

# eof
