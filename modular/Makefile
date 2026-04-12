CC      = gcc
CFLAGS  = -std=c11 -Wall -Wextra -O2

SRCS_COMMON = cusp.c p1.c manin.c sparse.c linalg.c hecke.c main.c
SRCS_TEST   = cusp.c p1.c manin.c sparse.c linalg.c hecke.c test_p1.c

# Default backend
BACKEND ?= builtin

ifeq ($(BACKEND),gmp)
  RAT_SRC   = rat_gmp.c
  RAT_DEF   = -DRAT_BACKEND_GMP
  LIBS      = -lgmp
  BIN       = modular_gmp
  TBIN      = test_p1_gmp
else ifeq ($(BACKEND),fp)
  RAT_SRC   = rat_fp.c
  RAT_DEF   = -DRAT_BACKEND_FP
  LIBS      =
  BIN       = modular_fp
  TBIN      = test_p1_fp
else
  RAT_SRC   = rat_builtin.c
  RAT_DEF   = -DRAT_BACKEND_BUILTIN
  LIBS      =
  BIN       = modular
  TBIN      = test_p1
endif

all: $(BIN)

$(BIN): $(SRCS_COMMON) $(RAT_SRC) rat.h
	$(CC) $(CFLAGS) $(RAT_DEF) -o $@ $(SRCS_COMMON) $(RAT_SRC) $(LIBS)

test: $(TBIN)
	./$(TBIN)

$(TBIN): $(SRCS_TEST) $(RAT_SRC) rat.h
	$(CC) $(CFLAGS) $(RAT_DEF) -o $@ $(SRCS_TEST) $(RAT_SRC) $(LIBS)

clean:
	rm -f modular modular_gmp modular_fp test_p1 test_p1_gmp test_p1_fp *.o

.PHONY: all test clean
