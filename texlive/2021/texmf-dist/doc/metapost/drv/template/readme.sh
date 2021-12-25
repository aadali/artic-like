# In order to compile the sources and to get a pdf file,
# use `make' or `./readme.sh' or:

# 1. run metapost THREE TIMES:
mpost template.mp
mpost template.mp
mpost template.mp

# 2. run latex (once):
latex template.tex

# 3. run dvips (once):
dvips template.dvi -o

# 4. run ps2pdf (once):
ps2pdf template.ps
