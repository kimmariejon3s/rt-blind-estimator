#### Makefile

CC=gcc

BIN_NAME=blind_rt_est
C_FILES=main.c

CFLAGS=-o $(BIN_NAME) $(C_FILES) -march=native
LIB_PATH=./libs
INC_PATH=./include
INCLUDE_FLAGS=-I$(INC_PATH)

LD_DIRS=-L$(LIB_PATH)
LDFLAGS=-lsndfile -lswresample -lavutil -lliquid -lnlopt -lm -lsvd

FULL_PATH=$(shell pwd)

all:
	@echo "Building application..."
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) $(LD_DIRS) $(LDFLAGS)
	@echo "Done!"

clean:
	@echo "Cleaning up dependencies..."
	rm -rf $(BIN_NAME) 
