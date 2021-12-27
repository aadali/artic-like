#!/bin/bash

# this is how to run in case if style file interfers with alphabetical sorting:
# iconv -f UTF-8 -t CP1251 rumakeindex-ex-x.idx | makeindex -i -s shipunov2i.ist | iconv -f CP1251 -t UTF-8 > rumakeindex-ex-x.ind
# however, without style file, everything is much simpler:
makeindex rumakeindex-ex-x.idx
