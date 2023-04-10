CC      = clang
CFLAGS  = -Ofast -march=native -flto=thin 
INLINE_FLAGS = #-finline-aggressive -finline-recursion=4
LOOP_FLAGS = -faggressive-loop-transform #-funroll-all-loops 
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
