
FP16 := #on
FP32 := #on
FP64 := on

SRCDIR := ./src
OBJDIR := ./obj
BINDIR := ./bin

CCHOME := /usr
CUDAHOME := /public/software/cuda-11.5
MPIHOME := /public/software/openmpi-4.1.1-cuda.10

CC := $(CCHOME)/bin/gcc -pipe
CXX := $(CCHOME)/bin/g++ -pipe
GC := $(CUDAHOME)/bin/nvcc -rdc=true -maxrregcount=127 -arch=sm_80 #-Xptxas=-v

GCC := $(CXX)

INCS := -I ./inc

LIBS := -L$(CUDAHOME)/lib64 -lcudart -lcublas
INCS += -I$(CUDAHOME)/include 

LIBS += -L$(MPIHOME)/lib -lmpi
INCS += -I$(MPIHOME)/include 

LIBS += -lm

CFLAGS := -c -O2 -std=c++11
LFLAGS := -O2

GCFLAGS := 
# GCFLAGS += -g -G

GCFLAGS += -x cu

vpath

vpath % $(SRCDIR)
vpath % $(OBJDIR)
vpath % $(BINDIR)

DFLAGS_LIST := FP16 FP32 FP64

DFLAGS := $(foreach flag,$(DFLAGS_LIST),$(if $($(flag)),-D$(flag)))

OBJS := main.o possion.o common_func.o iterative.o

OBJS := $(addprefix $(OBJDIR)/,$(OBJS))


$(BINDIR)/solver: $(OBJS)
	$(GCC) $(LFLAGS) $(LIBS) $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GCC) $(CFLAGS) $(DFLAGS) $(INCS)  $^ -o $@

clean:
	-rm $(OBJDIR)/* -rf
	-rm $(BINDIR)/solver

