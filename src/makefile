CFLAGS = -w -Wall -O3
BASIC_DEPS_H = bitmap.h profile_helper.h constants.hh shapes.hh kvec.hh options.hh
BASIC_DEPS_C = obstack.o bitmap.o profile_helper.o
PESTRIE_DEPS_H = pestrie.hh segtree.hh treap.hh histogram.hh matrix-ops.hh
PESTRIE_DEPS_C = segtree.o treap.o pes-common.o pes-self.o pes-dual.o matrix-ops.o
BITINDEX_DEPS_H = matrix-ops.hh bit-index.hh
BITINDEX_DEPS_C = matrix-ops.o bit-pt.o bit-se.o
LIB = #-L/usr/local/lib -ltcmalloc
CC = g++

ifndef INSTALL_DIR
	INSTALL_DIR=$(HOME)
endif


all: pesI bitI qtester


obstack.o: obstack.cc
	$(CC) obstack.cc $(CFLAGS) $(LIB) -c

bitmap.o: bitmap.cc bitmap.h obstack.o
	$(CC) bitmap.cc $(CFLAGS) $(LIB) -c

profile_helper.o: profile_helper.h profile_helper.cc
	$(CC) profile_helper.cc $(CFLAGS) $(LIB) -c

treap.o : treap.hh treap.cc $(BASIC_DEPS_H)
	$(CC) treap.cc $(CFLAGS) $(LIB) -c

segtree.o : segtree.hh segtree.cc $(BASIC_DEPS_H)
	$(CC) segtree.cc $(CFLAGS) $(LIB) -c

pes-common.o : pestrie.hh pes-common.cc $(BASIC_DEPS_H)
	$(CC) pes-common.cc $(CFLAGS) $(LIB) -c

pes-self.o : segtree.hh pes-self.cc $(BASIC_DEPS_H)
	$(CC) pes-self.cc $(CFLAGS) $(LIB) -c

pes-dual.o : segtree.hh pes-dual.cc $(BASIC_DEPS_H)
	$(CC) pes-dual.cc $(CFLAGS) $(LIB) -c

matrix-ops.o : matrix-ops.hh matrix-ops.cc
	$(CC) matrix-ops.cc $(CFLAGS) $(LIB) -c

bit-pt.o : bit-pt.cc bit-index.hh $(BASIC_DEPS_H)
	$(CC) bit-pt.cc $(CFLAGS) $(LIB) -c

bit-se.o : bit-se.cc bit-index.hh $(BASIC_DEPS_H)
	$(CC) bit-se.cc $(CFLAGS) $(LIB) -c

pes-querier.o : pes-querier.cc shapes.hh query.hh query-inl.hh options.hh $(BASIC_DEPS_H) $(BASIC_DEPS_C)
	$(CC) pes-querier.cc $(CFLAGS) $(LIB) -c

bit-querier.o : bit-querier.cc query.hh query-inl.hh options.hh $(BASIC_DEPS_H) $(BASIC_DEPS_C)
	$(CC) bit-querier.cc $(CFLAGS) $(LIB) -c

pesI: pes-indexer.cc  $(BASIC_DEPS_H) $(PESTRIE_DEPS_H) $(PESTRIE_DEPS_C) $(BASIC_DEPS_C)
	$(CC) pes-indexer.cc $(PESTRIE_DEPS_C) $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o pesI

bitI: bit-indexer.cc $(BASIC_DEPS_C) $(BASIC_DEPS_H) $(BITINDEX_DEPS_C) $(BITINDEX_DEPS_H)
	$(CC) bit-indexer.cc $(BASIC_DEPS_C) $(BITINDEX_DEPS_C) $(CFLAGS) $(LIB) -o bitI

qtester: qtester.cc pes-querier.o bit-querier.o matrix-ops.o query.hh options.hh $(BASIC_DEPS_H) $(BASIC_DEPS_C)
	$(CC) qtester.cc pes-querier.o bit-querier.o matrix-ops.o $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o qtester

formatter: matrix-ops.hh matrix-ops.cc formatter.cc
	$(CC) formatter.cc matrix-ops.o $(BASIC_DEPS_C) -o formatter

install:
	cp pesI bitI qtester $(INSTALL_DIR)/bin

clean:
	rm -f *.o pesI bitI qtester formatter

