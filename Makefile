CC      = gcc
CFLAGS  = -Ofast -march=core-avx2 -ftree-vectorize -fopt-info-vec -fopt-info-vec-missed -fopenmp
WFLAGS  = -std=c11 -Wall -Wextra -Werror -g
LDFLAGS = -lm

TARGETS = tiny_md viz
SOURCES = $(shell echo *.c)
OBJECTS = core.o wtime.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(WFLAGS) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all                      
