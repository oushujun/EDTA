CC = gcc
CXX = g++

# For program run
# CPPFLAGS = -fomit-frame-pointer -W -Wall -Winline -DNDBUG -O3
CPPFLAGS = -fomit-frame-pointer -DNDBUG  -O3 
CFLAGS = -fomit-frame-pointer -DNDBUG  -O3 

# For debug
# CPPFLAGS = -g -W -Wall -Winline

ALL_OPS = regex.o Timing.o BooleanString.o LinearSuffixSort.o lcp.o \
          stdaln.o seq.o PairsFilter.o PBS.o CallPsScan.o \
          stdaln_interface.o

ALL:$(ALL_OPS) ltr_finder psearch

regex.o:regex.c
	$(CC) -c regex.c  $(CFLAGS)
	
stdaln_interface.o:stdaln_interface.cpp
	$(CXX) -c stdaln_interface.cpp $(CPPFLAGS)

CallPsScan.o:CallPsScan.cpp
	$(CXX) -c CallPsScan.cpp $(CPPFLAGS)
	
PBS.o: PBS.cpp
	$(CXX) -c PBS.cpp $(CPPFLAGS)

PairsFilter.o: PairsFilter.cpp
	$(CXX) -c PairsFilter.cpp $(CPPFLAGS)

lcp.o:lcp.c
	$(CXX) -c lcp.c $(CPPFLAGS)

BooleanString.o:BooleanString.C BooleanString.h
	$(CXX) -c BooleanString.C $(CPPFLAGS)

Timing.o:Timing.C
	$(CXX) -c Timing.C $(CPPFLAGS)

LinearSuffixSort.o:LinearSuffixSort.cpp
	$(CXX) -c LinearSuffixSort.cpp $(CPPFLAGS)

stdaln.o:stdaln.c
	$(CXX) -c stdaln.c $(CPPFLAGS)
seq.o:seq.c
	$(CXX) -c seq.c $(CPPFLAGS)

ltr_finder: $(ALL_OPS) main.cpp
	$(CXX) -o ltr_finder main.cpp PairsFilter.o stdaln.o seq.o LinearSuffixSort.o BooleanString.o Timing.o lcp.o PBS.o CallPsScan.o stdaln_interface.o regex.o $(CPPFLAGS)

psearch: $(ALL_OPS) psearch.cpp
	$(CXX) -o psearch psearch.cpp PBS.o seq.o stdaln_interface.o stdaln.o

clean:
	rm -f *.o





