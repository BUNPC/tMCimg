SRCS = tMCimgLOT.c
OBJS = $(SRCS:%.c=%.o)

default: opt

tMCimgLOT: $(OBJS) makefile
	gcc  -pedantic -Wall $(LINK_OPTIONS) -o tMCimgLOT $(OBJS) -lm

%.o: %.c makefile
	gcc  -pedantic -Wall $(COMPILE_OPTIONS) -c  $< -o $@

# debug: COMPILE_OPTIONS = -g -DDEBUG -DSINGLE_PREC -DNOT_POINTSRC
debug: COMPILE_OPTIONS = -g -DDEBUG -DSINGLE_PREC -DREFLECT
debug: LINK_OPTIONS = -g 
debug: tMCimgLOT

# opt:    COMPILE_OPTIONS = -O3 -DSINGLE_PREC -DNOT_POINTSRC
opt:    COMPILE_OPTIONS = -O3 -DSINGLE_PREC -DREFLECT
double: COMPILE_OPTIONS = -O3 

double opt: LINK_OPTIONS = -O3
double opt: tMCimgLOT



