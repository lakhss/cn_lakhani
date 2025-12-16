#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include/openblas
OPTCLOCAL=-fPIC -march=native
