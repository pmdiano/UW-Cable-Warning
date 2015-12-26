CC = g++

CFLAGS = `pkg-config --cflags opencv` -Wall --std=c++0x

LDFLAGS = -lfftw3 -lm `pkg-config --libs opencv`

SOURCES = src/cable.cc src/uwcw.cc src/particles.cc

OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = uwcw.out

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJECTS)

rebuild: clean all
