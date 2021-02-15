## Makefile

PGM	= IDAstar

all: $(PGM).c
	   gcc -Wall  -o $(PGM) $(PGM).c -I.


# -Wl,--stack,1677216
# -fconserve-stack