SHELL=/bin/zsh
CC=gcc
CFLAGS=-g -c -fdiagnostics-color=always -std=c99 -Wall -MMD -MP
LIBFLAGS=-lm
ELF = lab2

.PHONY: clean all rebuild

all : $(ELF)

$(ELF) : $(patsubst %.c, %.o, $(shell print **/*.c))
	$(CC) $(LIBFLAGS) $^ -o $@

%.o : %.c
	$(CC) $(CFLAGS) $< -o $@

clean: clean_specific
	rm -f **/*.o(N) **/*.d(N) $(ELF)

clean_specific:


rebuild: clean all

-include $(patsubst %.c, %.d, $(shell print **/*.c))
