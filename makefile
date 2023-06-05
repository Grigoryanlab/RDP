# flags
CC := g++
CPP_FLAGS := -std=c++11 -fPIC
MST_INC_DIR :=
MST_OBJS_DIR :=
MST_OJBS_LIST := mstcondeg.o mstoptions.o mstrotlib.o msttransforms.o msttypes.o mstfasst.o mstsequence.o
MST_OBJS := $(foreach dep,$(MST_OJBS_LIST),$(MST_OBJS_DIR)/$(dep))

all : mstObjsDir.txt dockingDistribution.o dockingDistribution

mstObjsDir.txt :
	@echo $(MST_OBJS_DIR) > mstObjsDir.txt

dockingDistribution.o : src/dockingDistribution.cpp
	@mkdir -p objs 
	$(CC) $(CPP_FLAGS) -c -o objs/dockingDistribution.o -I$(MST_INC_DIR) src/dockingDistribution.cpp

dockingDistribution : objs/dockingDistribution.o src/dockingDistribution.cpp
	@mkdir -p bin 
	@mkdir -p matchedBackbones 
	@mkdir -p matchedBackbones/P1
	@mkdir -p matchedBackbones/P2
	$(CC) $(CPP_FLAGS) -o bin/dockingDistribution  -Iinclude $(MST_OBJS) objs/dockingDistribution.o
