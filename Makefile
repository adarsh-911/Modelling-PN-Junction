# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O2 -Iinc

# Source files
SRCS = main.c \
       src/equillibrium.c \
       src/init.c \
       src/poisson_solver.c \
			 src/steady_state/non_linear_newton.c \
			 src/steady_state/steady_state.c

# Object files
OBJS = $(SRCS:.c=.o)

# Executable name
TARGET = main

# Default rule
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

# Compile .c to .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

