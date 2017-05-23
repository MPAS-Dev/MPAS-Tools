.SUFFIXES: .c .o


OBJS = edgelist.o geometry.o heap.o memory.o output.o voronoi.o voronoi_main.o

all: $(OBJS)
	ar -ru libfortune.a $(OBJS)

clean:
	rm -f *.o libfortune.a

.c.o:
	$(CC) $(CFLAGS) -c $< -I../
