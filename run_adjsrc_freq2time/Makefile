CC = gcc
CFLAGS = -g -Wall


BIN = .
INC = -I.
LIB =  -lm
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)



main: clean $(OBJ) 
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)


%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean :
	rm -f $(BIN)/main *~ *.o basis_function
