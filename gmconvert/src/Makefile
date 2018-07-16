CC	      = gcc

CFLAGS	      = -O3 -g -Wall

DEST	      = .

EXTHDRS	      =

HDRS	      = \

INSTALL	      = /etc/install

LD	      = gcc

LDFLAGS	      = -O3 -g -Wall -lm

LIBS	      = 

MAKEFILE      = Makefile

OBJS	      = gmconvert.o\
	        Radius.o\
	        Voxel.o\
	        MCubeFunc.o\
	        MCubeIO.o\
	        MCubeVolSur.o\
	        GridMap.o\
	        PointEM.o\
	        GmmAtomEM.o\
	        AtomKmean.o\
	        3DmapKmean.o\
	        3DmapEM.o\
	        Gmm3DmapEM.o\
	        MapCCP4.o\
	        GaussIO.o\
	        BasicIO.o\
	        Atom2Vox.o\
	        TabExpMinX.o\
	        qRMS.o\
	        Matrix3D.o\
	        HeapSort.o\
	        Jacobi3.o\
    	        PdbIO.o \
    	        io_mmCIF.o \
    	        PDB_from_mmCIF.o \
    	        ATOMs_gmconvert_from_PDBmmCIF.o \
    	        Ellipsoid.o

PRINT	      = pr

PROGRAM       = ../gmconvert 

SHELL	      = /bin/sh

SRCS	      = gmconvert.c \
	        Radius.c \
	        Voxel.c \
	        MCubeFunc.c \
	        MCubeIO.c \
	        MCubeVolSur.c \
	        GridMap.c \
	        PointEM.c \
	        GmmAtomEM.c \
	        AtomKmean.c \
	        3DmapKmean.c \
	        3DmapEM.c \
	        Gmm3DmapEM.c \
	        MapCCP4.c \
	        GaussIO.c \
	        BasicIO.c \
	        Atom2Vox.c \
	        TabExpMinX.c \
	        qRMS.c \
	        Matrix3D.c \
	        HeapSort.c \
	        Jacobi3.c \
    	        PdbIO.c \
    	        io_mmCIF.c \
    	        PDB_from_mmCIF.c \
    	        ATOMs_gmconvert_from_PDBmmCIF.c \
    	        Ellipsoid.c

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(OBJS) $(LDFLAGS) $(LIBS)  -o $(PROGRAM)
		@echo "done"

