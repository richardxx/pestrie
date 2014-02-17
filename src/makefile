CFLAGS = -w -Wall -O2
BASIC_DEPS_H = bitmap.h profile_helper.h options.h 
BASIC_DEPS_C = obstack.o bitmap.o profile_helper.o
PESTRIE_DEPS_H = pestrie.hh segtree.hh treap.hh kvec.hh histogram.hh matrix_ops.hh
PESTRIE_DEPS_C = segtree.o treap.o pes-common.o pes-self.o pes-dual.o matrix_ops.o
BITINDEX_DEPS_H = matrix_ops.hh bit_index.hh
BITINDEX_DEPS_C = matrix_ops.o bit_pt.o bit_se.o
LIB = #-ltcmalloc
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

pes-common.o : pestrie.hh pes-common.cc
	$(CC) pes-common.cc $(CFLAGS) $(LIB) -c

pes-self.o : segtree.hh pes-self.cc
	$(CC) pes-self.cc $(CFLAGS) $(LIB) -c

pes-dual.o : segtree.hh pes-dual.cc
	$(CC) pes-dual.cc $(CFLAGS) $(LIB) -c

matrix_ops.o : matrix_ops.hh matrix_ops.cc
	$(CC) matrix_ops.cc $(CFLAGS) $(LIB) -c

bit_pt.o : bit_pt.cc bit_index.hh
	$(CC) bit_pt.cc $(CFLAGS) $(LIB) -c

bit_se.o : bit_se.cc bit_index.hh
	$(CC) bit_se.cc $(CFLAGS) $(LIB) -c

pesI: pes-indexer.cc  $(BASIC_DEPS_H) $(PESTRIE_DEPS_H) $(PESTRIE_DEPS_C) $(BASIC_DEPS_C)
	$(CC) pes-indexer.cc $(PESTRIE_DEPS_C) $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o pesI

bitI: bit-indexer.cc $(BASIC_DEPS_C) $(BASIC_DEPS_H) $(BITINDEX_DEPS_C) $(BITINDEX_DEPS_H)
	$(CC) bit-indexer.cc $(BASIC_DEPS_C) $(BITINDEX_DEPS_C) $(CFLAGS) $(LIB) -o bitI

pesQ: pes-querier.cc treap.hh kvec.hh histogram.hh query.h $(BASIC_DEPS_H) $(BASIC_DEPS_C)
	$(CC) pes-querier.cc $(BASIC_DEPS_C) $(CFLAGS) $(LIB) -o pesQ

bitQ: bit-querier.cc query.h $(BASIC_DEPS_C) $(BASIC_DEPS_H) matrix_ops.o
	$(CC) bit-querier.cc $(BITINDEX_DEPS_C) $(BASIC_DEPS_C) $(BASIC_DEPS_H) $(CFLAGS) $(LIB) -o bitQ

formatter: matrix_ops.hh matrix_ops.cc formatter.cc
	$(CC) formatter.cc matrix_ops.o $(BASIC_DEPS_C) -o formatter

install:
	cp pesI bitI pesQ bitQ $(INSTALL_DIR)/bin

clean:
	rm -f *.o pesI bitI pesQ bitQ formatter

