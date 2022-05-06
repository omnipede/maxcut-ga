OUT = maxcut
CC = gcc
SDIR = .
CFLAGS = -Wall -c

_OBJS = main.o graph.o util.o
OBJS = $(patsubst %,$(SDIR)/%,$(_OBJS))

$(SDIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c $(INC) -o $@ $< $(CFLAGS)

all: $(OBJS)
	$(CC) $(OBJS) -o $(OUT)

run:
	./$(OUT) ./maxcut.in ./maxcut.out

.PHONY: clean

clean:
	rm -f *.o
	rm $(OUT)