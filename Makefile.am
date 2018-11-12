#INCLUDES = syspara.h

bin_PROGRAMS = KRT
#KRT_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
KRT_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -vec-report2 -pg

KRT_INCLUDE = -I./include
KRT_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

KRT_SOURCES = ./src/krt.c ./src/val_constants.c ./src/mem.c\
	./src/win_out.c ./src/input_para.c ./src/system.c ./src/currents.c\
	./src/eventloop.c ./src/exptable.c\
	./src/eular.c ./src/runge.c ./src/syspara.h

