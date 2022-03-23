all: fqzcomp5

INCLUDES=-Ihtscodecs
LIBS=-Lhtscodecs/htscodecs/.libs -lhtscodecs
LDFLAGS=
CFLAGS=-g -O3
#CFLAGS=-g

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

OBJ=fqzcomp5.o lzp16e.o
fqzcomp5: htscodecs $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $@ $(LIBS) -pthread -lm -lbz2

htscodecs:
	cd htscodecs; $(MAKE)

clean:
	-rm fqzcomp5 *.o
