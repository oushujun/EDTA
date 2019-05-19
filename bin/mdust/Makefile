# Useful directories

# Directories to search for header files
SEARCHDIRS := -I. -I- 


#SYSTYPE :=     $(shell uname)

CC      := gcc
CFLAGS  = -O2 -Wall ${SEARCHDIRS} -D_REENTRANT

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := gcc
LDFLAGS    =
LOADLIBES := 

.PHONY : all
all:    mdust

mdust:  ./mdust.o ./dust.o ./fastafile.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBES}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core *.o 

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core *.o 


