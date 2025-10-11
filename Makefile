CC = gcc
CFLAGS = -O2 -Iinc

SRCS = main.c \
       src/equillibrium.c \
       src/init.c \
       src/poisson_solver.c \
			 src/steady_state/non_linear_newton.c \
			 src/steady_state/steady_state.c

OBJS = $(SRCS:.c=.o)

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

