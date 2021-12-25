# In order to compile the sources and to get a pdf file,
# use `make' or `./readme.sh' or:

# 1. run metapost TWICE:
mpost coq-sample.mp
mpost coq-sample.mp

# 2. run latex (twice):
latex coq-sample.tex
latex coq-sample.tex

# 3. run dvips (once):
dvips coq-sample.dvi -o

# 4. run ps2pdf (once):
ps2pdf coq-sample.ps
