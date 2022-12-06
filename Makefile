CC	  = gcc
AR	  = ar
CFLAGS	  = $(INCLUDES) $(DEFINES) $(DEBUG) $(OPTIMIZE)

OBJS  = fft1d.o fft2d.o fftmem.o

LIBS	  = -lm

TARGET	  =  libfftlib.a

DEBUG	  = -Wall #-g -Wall
DEFINES	  = 
INCLUDES  = 
OPTIMIZE  =


###-----------------------------------------------

all: $(TARGET)

clean:
	rm -f *.o
	rm -f core
	rm -f *~ 
	rm $(TARGET)

$(TARGET) : $(OBJS)
	$(AR) rs $(TARGET) $(OBJS)


##--- Compile
.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -c $<

fft1d.o : fft1d.c 
fft2d.o : fft2d.c
fftmem.o : fftmem.c