CC = gcc

CC_FLAGS = -Wall -Wextra -g -O0
LDD_FLAGS = -lm

OUTPUT = main.out
SRC = main.c

.PHONY = run

all: $(OUTPUT)

run: $(OUTPUT)
	@./$(OUTPUT)

$(OUTPUT): $(SRC)
	$(CC) $(LDD_FLAGS) $^ -o $@

