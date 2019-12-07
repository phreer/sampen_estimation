WARNS += -Wall -Wextra
STD = -std=c++11
DEBUG_FLAG += -g -DDEBUG
RELEASE_FLAG += -O3

OBJS = random_sampler.o

LIBS_DIR += -L$(HOME)/local/lib
LIBS += -lgsl
INCLUDE_PATH += -I$(HOME)/local/include
CXXFLAGS += $(WARNS) $(STD) $(LIBS_DIR) $(LIBS) ${INCLUDE_PATH}

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS += $(DEBUG_FLAG)
else
	CXXFLAGS += $(RELEASE_FLAG)
endif

random_sampler.o: random_sampler.cpp random_sampler.h tensor.h
	$(CXX) ${CXXFLAGS} -c $<

sampen_range_tree: sampen_range_tree.cpp  RangeTree2.h $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJS)
	echo "DEBUG mode = $(DEBUG)"

clean:
	rm sampen_range_tree *.o