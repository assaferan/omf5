VPATH=src
CC=clang
CFLAGS=-I${GMP_PATH}/include -I${ANTIC_PATH}/include -Iinclude
LDFLAGS=-L${GMP_PATH}/lib -L${ANTIC_PATH}/lib -lgmp -lflint -lantic -lfunctions

omf5 : omf5.o arith.o fq_nmod_mat.o fq_nmod_mpoly.o fq_nmod_mpoly_mat.o fq_nmod_quad.o genus.o hash.o hecke.o jordan.o mass.o matrix_tools.o nbr_data.o neighbor.o nf_mat.o orthogonalize.o pivot_data.o qf_inv.o tests.o
	$(CC) $(LDFLAGS) $^ -o $@
