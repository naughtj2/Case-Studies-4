CC = gcc
CFLAGS = -Wall -O2 

TARGETS = Question1 Question2 Question3 Question3b

all: $(TARGETS)

Question1: Question1.c
	$(CC) $(CFLAGS) -o $@ $< -lm

Question2: Question2.c
	$(CC) $(CFLAGS) -o $@ $< -lm

Question3: Question3.c
	$(CC) $(CFLAGS) -o $@ $< -lm

Question3b: Question3b.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f $(TARGETS)

