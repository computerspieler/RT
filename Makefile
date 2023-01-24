COMMON_SRCS=$(wildcard common/*.c)
SRCS=$(wildcard src/*.c)
FULL_SRCS=$(SRCS) $(COMMON_SRCS)
DEPS=bin/deps/kernel/main.cl.d $(patsubst %,bin/deps/%.d,$(FULL_SRCS))
OBJS=$(patsubst %,bin/objs/%.o,$(FULL_SRCS))

CC=gcc
CPP=gcc -E
LD=gcc
CCFLAGS=-Wall -Wextra -Isrc -Ibin/headers -Iinclude -Icommon -c \
	-Wno-incompatible-pointer-types
CPPFLAGS=-Iinclude -P -nostdinc -D IN_OPENCL -D PRINT_DEBUG -fdirectives-only
LDFLAGS=-lm `sdl-config --cflags --libs` -lOpenCL

all: debug

clean:
	rm -rf bin

-include $(DEPS)

debug: bin/debug
release: bin/release

bin/debug: CCFLAGS := $(CCFLAGS) -g -ggdb
bin/debug: LDFLAGS := $(LDFLAGS) -g -ggdb
bin/debug: bin/headers/kernel/main.cl.h $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

bin/release: CCFLAGS := $(CCFLAGS) -O3
bin/release: LDFLAGS := $(LDFLAGS) -O3
bin/release: bin/headers/kernel/main.cl.h $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

bin/headers/%.cl.h: %.cl
	@mkdir -p $(dir $@)
	
	echo "char $(subst .,_,$(subst /,_,$<))_src[] = " > $@
	cat $< $(COMMON_SRCS) | $(CPP) $(CPPFLAGS) - | sed 's/\\/\\\\/g' | sed 's/\"/\\\"/g' | sed -e 's/\(.*\)/"\1\\n"/' >> $@
	echo ";" >> $@

bin/deps/%.cl.d: %.cl
	@mkdir -p $(dir $@)
	$(CPP) $(CPPFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%.d,bin/headers/%.h,$@)

bin/objs/%.c.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -o $@ $<

bin/deps/%.c.d: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%.d,bin/objs/%.o,$@)

