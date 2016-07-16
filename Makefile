#### Makefile

CC=gcc

CFLAGS=-o blind_rt_est main.c -march=native
EXT_LIB_PATH=./lib_src
INCLUDE_FLAGS=-I./ -I$(EXT_LIB_PATH)/SVDLIBC

LD_DIRS=-L$(EXT_LIB_PATH)/libsndfile-1.0.27/src/.libs/libsndfile.so
LD_DIRS+=-L$(EXT_LIB_PATH)/SVDLIBC

LDFLAGS=-lsndfile -lswresample -lavutil -lliquid -lnlopt -lm -lsvd

FULL_PATH=$(shell pwd)
SVDLIBC_PATH=$(EXT_LIB_PATH)/SVDLIBC
NLOPT_PATH=$(EXT_LIB_PATH)/nlopt-2.4.2
FFMPEG_PATH=$(EXT_LIB_PATH)/ffmpeg-3.1.1
LIBSND_PATH=$(EXT_LIB_PATH)/libsndfile-1.0.27
LIQUID_PATH=$(EXT_LIB_PATH)/liquid-dsp-1.2.0

all: extract build_deps compile

compile:
	@echo "Building application..."
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) $(LD_DIRS) $(LDFLAGS)
	@echo "Done!"

build_deps: extract
	@echo "Building dependencies..."
	cd $(SVDLIBC_PATH); make
	cd $(NLOPT_PATH); ./configure --prefix="$(FULL_PATH)/libs"; make; make install
	cd $(FFMPEG_PATH); ./configure --prefix="$(FULL_PATH)/libs" --bindir="$(FULL_PATH)/libs"; make; make install 
	cd $(LIBSND_PATH); ./configure; make
	cd $(LIQUID_PATH); ./reconf; ./configure --prefix="$(FULL_PATH)/libs"; make; make install

extract:
	cd $(EXT_LIB_PATH); make

clean:
	@echo "Cleaning up dependencies..."
	cd $(SVDLIBC_PATH); make clean
	cd $(NLOPT_PATH); make clean
	cd $(FFMPEG_PATH); make clean
	cd $(LIBSND_PATH); make clean
	cd $(LIQUID_PATH); make clean
	cd $(EXT_LIB_PATH); make clean
