# the compiler: gcc for C program, define as g++ for C++
CXX = g++
RM = rm -f

# compiler flags:
CFLAGS  = -O3 -std=c++0x
LDFLAGS =

# the build target executable:
SOURCES = MostCommonK-mers.cpp
SOURCE_DIR = src
IN_SOURCES = $(addprefix $(SOURCE_DIR)/, $(SOURCES))

OBJECTS = $(subst .cpp,.o, $(SOURCES))
EXEC = MostCommonK-mers

all: $(OBJECTS)
	$(CXX) $(CFLAGS) -o $(EXEC) $(IN_SOURCES)

$(OBJECTS):
	$(CXX) $(CFLAGS) -c $(IN_SOURCES)

clean:
	$(RM) $(OBJECTS)
