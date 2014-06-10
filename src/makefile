CFLAGS = -w -Wall -o3
BASIC_DEPS_H = bitmap.h profile_helper.h options.h 
BASIC_DEPS_C = obstack.o bitmap.o profile_helper.o
PESTRIE_DEPS_H = pestrie.hh segtree.hh treap.hh kvec.hh histogram.hh matrix-ops.hh
PESTRIE_DEPS_C = segtree.o treap.o pes-common-exp.o pes-self.o pes-dual.o matrix-ops.o
BITINDEX_DEPS_H = matrix-ops.hh bit-index.hh
BITINDEX_DEPS_C = matrix-ops.o bit-pt.o bit-se.o
LIB = #-L/usr/local/lib -ltcmalloc
CC = g++

ifndef INSTALL_DIR
	INSTALL_DIR=$(HOME)
endif


all: pesI bitI pesQ bitQ


obstack.o: obstack.cc
	$(CC) obstack.cc $(CFLAGS) $(LIB) -c

bitmap.o: bitmap.cc bitmap.h obstack.o
	$(CC) bitmap.cc $(CFLAGS) $(LIB) -c

profile_helper.o: profile_helper.h profile_helper.cc
	$(CC) profile_helper.cc $(CFLAGS) $(LIB) -c

treap.o : treap.hh treap.cc
	$(CC) treap.cc $(CFLAGS) $(LIB) -c

segtree.o : segtree.hh segtree.cc
	$(CC) segtree.cc $(CFLAGS) $(LIB) -c

pes-common-exp.o : pestrie.hh pes-common-exp.cc
	$(CC) pes-common-exp.cc $(CFLAGS) $(LIB) -c

pes-self.o : segtree.hh pes-self.cc
	$(CC) pes-self.cc $(CFLAGS) $(LIB) -c

pes-dual.o : segtree.hh pes-dual.cc
	$(CC) pes-dual.cc $(CFLAGS) $(LIB) -c

matrix-ops.o : matrix-ops.hh matrix-ops.cc
	$(CC) matrix-ops.cc $(CFLAGS) $(LIB) -c

bit-pt.o : bit-pt.cc bit-index.hh
	$(CC) bit-pt.cc $(CFLAGS) $(LIB) -c

bit-se.o : bit-se.cc bit-index.hh
	$(CC) bit-se.cc $(CFLAGS) $(LIB) -c

pesI: pes-indexer.cc  $(BASIC_DEPS_H) $(PESTRIE_DEPS_H) $(PESTRIE_DEPS_C) $(BASIC_DEPS_C)
	$(CC) pes-indexer.cc $(PESTRIE_DEPS_C) $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o pesI

bitI: bit-indexer.cc $(BASIC_DEPS_C) $(BASIC_DEPS_H) $(BITINDEX_DEPS_C) $(BITINDEX_DEPS_H)
	$(CC) bit-indexer.cc $(BASIC_DEPS_C) $(BITINDEX_DEPS_C) $(CFLAGS) $(LIB) -o bitI

pesQ: pes-querier-exp.cc treap.hh kvec.hh histogram.hh query.h pes-querier.hh $(BASIC_DEPS_H) $(BASIC_DEPS_C)
	$(CC) pes-querier-exp.cc $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o pesQ

bitQ: bit-querier.cc query.h $(BASIC_DEPS_C) $(BASIC_DEPS_H) matrix-ops.o
	$(CC) bit-querier.cc $(BITINDEX_DEPS_C) $(BASIC_DEPS_C) $(BASIC_DEPS_H) $(CFLAGS) $(LIB) -o bitQ

formatter: matrix-ops.hh matrix-ops.cc formatter.cc
	$(CC) formatter.cc matrix-ops.o $(BASIC_DEPS_C) -o formatter

install:
	cp pesI bitI pesQ bitQ $(INSTALL_DIR)/bin

clean:
	rm -f *.o pesI bitI pesQ bitQ formatter

