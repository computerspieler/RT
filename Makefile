SRCS_DIRS=$(sort $(dir $(wildcard src/*/))) src/ common/

SRCS=$(foreach DIR,$(SRCS_DIRS),$(wildcard $(DIR)*.c))
DEPS=$(patsubst %,bin/deps/%.d,$(SRCS))
OBJS=$(patsubst %,bin/objs/%.o,$(SRCS))

CC=clang
LD=clang
CCFLAGS=-Wall -Wextra -Isrc -Iinclude -Icommon -c \
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

bin/objs/%.c.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -o $@ $<

bin/deps/%.c.d: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%,bin/objs/%,$@)

