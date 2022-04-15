SRCS_DIRS=$(sort $(dir $(wildcard src/*/))) src/

SRCS=$(foreach DIR,$(SRCS_DIRS),$(wildcard $(DIR)*.c))
DEPS=$(patsubst src/%,bin/deps/%.d,$(SRCS))
OBJS=$(patsubst src/%,bin/objs/%.o,$(SRCS))

CC=gcc
LD=gcc
CCFLAGS=-Wall -Wextra $(patsubst %/,-I%,$(SRCS_DIRS)) -c \
	-Wno-incompatible-pointer-types
LDFLAGS=-lm `sdl-config --cflags --libs` -lOpenCL

all: debug

clean:
	rm -rf bin

-include $(DEPS)

debug: bin/debug
release: bin/release

bin/debug: CCFLAGS := $(CCFLAGS) -g -ggdb
bin/debug: LDFLAGS := $(LDFLAGS) -g -ggdb
bin/debug: $(OBJS)
	$(LD) -o $@ $^ $(LDFLAGS)

bin/release: CCFLAGS := $(CCFLAGS) -O2
bin/release: LDFLAGS := $(LDFLAGS) -O2
bin/release: $(OBJS)
	$(LD) -o $@ $^ $(LDFLAGS)

bin/objs/%.c.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -o $@ $<

bin/deps/%.c.d: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%,bin/objs/%,$@)

