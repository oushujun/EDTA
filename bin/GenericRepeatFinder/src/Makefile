# Define required macros here
CFLAGS = -Wall -std=c++11 -O3 -fopenmp
CC = g++
INCLUDES = -Igrf-main -Igrf-intersperse -Igrf-mite-cluster

all: grf-main grf-alignment grf-alignment2 grf-dbn grf-filter grf-intersperse grf-mite-cluster grf-nest

# grf-main
OBJS_GRF_MAIN = grf-main/main.o grf-main/DetectDR.o grf-main/DetectIR.o grf-main/DetectMITE.o grf-main/functions.o

grf-main: ${OBJS_GRF_MAIN}
	${CC} ${CFLAGS} ${OBJS_GRF_MAIN} -o $@/$@

# grf-alignment
OBJS_GRF_ALIGNMENT = grf-alignment/main.o

grf-alignment: ${OBJS_GRF_ALIGNMENT}
	${CC} ${CFLAGS} ${OBJS_GRF_ALIGNMENT} -o $@/$@

# grf-alignment2
OBJS_GRF_ALIGNMENT2 = grf-alignment2/main.o

grf-alignment2: ${OBJS_GRF_ALIGNMENT2}
	${CC} ${CFLAGS} ${OBJS_GRF_ALIGNMENT2} -o $@/$@

# grf-dbn
OBJS_GRF_DBN = grf-dbn/main.o

grf-dbn: ${OBJS_GRF_DBN}
	${CC} ${CFLAGS} ${OBJS_GRF_DBN} -o $@/$@

# grf-filter
OBJS_GRF_FILTER = grf-filter/main.o

grf-filter: ${OBJS_GRF_FILTER}
	${CC} ${CFLAGS} ${OBJS_GRF_FILTER} -o $@/$@

# grf-intersperse
OBJS_GRF_INTERSPERSE = grf-intersperse/main.o grf-intersperse/DetectIntersperse.o

grf-intersperse: ${OBJS_GRF_INTERSPERSE}
	${CC} ${CFLAGS} ${OBJS_GRF_INTERSPERSE} -o $@/$@

# grf-mite-cluster
OBJS_GRF_MITE_CLUESTER = grf-mite-cluster/main.o grf-mite-cluster/FamilyCluster.o

grf-mite-cluster: ${OBJS_GRF_MITE_CLUESTER}
	${CC} ${CFLAGS} ${OBJS_GRF_MITE_CLUESTER} -o $@/$@

# grf-nest
OBJS_GRF_NEST = grf-nest/main.o

grf-nest: ${OBJS_GRF_NEST}
	${CC} ${CFLAGS} ${OBJS_GRF_NEST} -o $@/$@

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $< -o $@

clean:
	rm -f grf-main/*.o
	rm -f grf-alignment/*.o
	rm -f grf-alignment2/*.o
	rm -f grf-dbn/*.o
	rm -f grf-filter/*.o
	rm -f grf-intersperse/*.o
	rm -f grf-mite-cluster/*.o
	rm -f grf-nest/*.o


