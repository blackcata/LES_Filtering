F90=ifort
FCFLAGS=-O3 -qopenmp#-I${MKL_HOME}/include/intel64/ilp64 -L${MKL_HOME}/lib -lmkl_intel_lp64   \
										 -lmkl_intel_thread -lmkl_core -lpthread

TARGET= LES_Filtering
OBJECT= LES_Filtering_module.o LES_Filtering_main.o LES_Filtering_setup.o \
				LES_Filtering_read.o LES_Filtering_filter.o LES_Filtering_second_filter.o \
				LES_Filtering_output.o LES_Filtering_eig33.o Vortical_Structure.o

all : $(TARGET)
$(TARGET) : $(OBJECT)
	$(F90) $(FCFLAGS) -o $@ $^

.SUFFIXES. : .o .f90

%.o : %.f90
	$(F90) $(FCFLAGS) -c $<

clean :
	rm -f *.o
	rm LES_Filtering les_filtering_module.mod
