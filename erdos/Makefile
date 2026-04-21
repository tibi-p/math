CC      = gcc
CFLAGS  = -O3 -march=native -std=c11 -Wall -Wextra

.PHONY: all clean run315 run945

all: dfs

dfs: dfs.c
	$(CC) $(CFLAGS) -o dfs dfs.c

run315: dfs
	./dfs 315

run945: dfs
	./dfs 945

clean:
	rm -f dfs
