CPP=g++
CPPFLAGS = -g -O2 -Wall -lboost_program_options
OBJ=readLengths.o

all: $(OBJ)
	$(CPP) $(CPPFLAGS) -o readLengths $(OBJ)

clean:
	rm -rf $(OBJ) readLengths