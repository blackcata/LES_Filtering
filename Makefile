F90=ifort
FCFLAGS=-O2

TARGET= LES_Filtering
OBJECT= LES_Filtering_module.o LES_Filtering_main.o LES_Filtering_setup.o \
				LES_Filtering_read.o

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean :
	rm -f *.o
