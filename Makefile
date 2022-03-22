all: fqzcomp5

INCLUDES=-Ihtscodecs
LIBS=-Lhtscodecs/htscodecs/.libs -lhtscodecs
LDFLAGS=
CFLAGS=-g -O3
#CFLAGS=-g

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

fqzcomp5: htscodecs fqzcomp5.o
	$(CC) $(LDFLAGS) fqzcomp5.o -o $@ $(LIBS) -pthread -lm -lbz2

htscodecs:
	cd htscodecs; $(MAKE)

clean:
	-rm fqzcomp5 *.o
