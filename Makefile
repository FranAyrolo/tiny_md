CC      = gcc
CFLAGS  = -O2 -march=native -ftree-vectorize -fopt-info-vec -fopt-info-vec-missed
WFLAGS  = -std=c11 -Wall -Wextra -Werror
LDFLAGS = -lm

TARGETS = tiny_md #viz
SOURCES = $(shell echo *.c)
OBJECTS = core.o wtime.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(WFLAGS) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all                      
