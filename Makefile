CC = g++

CFLAGS = -I/usr/local/opt/opencv3/include -Wall --std=c++0x -g

LDFLAGS = -lfftw3 -lm -L/usr/local/opt/opencv3/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs

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
