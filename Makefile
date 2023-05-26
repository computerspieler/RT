COMMON_SRCS=$(wildcard common/*.c)
SRCS=$(wildcard src/*.c)
FULL_SRCS=$(SRCS) $(COMMON_SRCS)
DEPS=$(patsubst %,bin/deps/%.d,$(FULL_SRCS))
OBJS=$(patsubst %,bin/objs/%.o,$(FULL_SRCS))
COMMON_OBJS=$(patsubst %,bin/objs/%.o,$(COMMON_SRCS))

PREFIX=
CC=$(PREFIX)gcc
CPP=$(PREFIX)gcc -E
LD=$(PREFIX)gcc
CCFLAGS=-Wall -Wextra -Isrc -Ibin/headers -Iinclude -Icommon -c -Wno-unused-parameter
CPPFLAGS=-Iinclude -P -nostdinc -D IN_OPENCL -fdirectives-only
LDFLAGS=-lm -lOpenCL -lGL -lglut -lGLU

all: debug

clean:
	rm -rf bin

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

debug: bin/debug
release: bin/release

bin/debug: CCFLAGS  := $(CCFLAGS)  -g -ggdb
bin/debug: LDFLAGS  := $(LDFLAGS)  -g -ggdb
bin/debug: CPPFLAGS := $(CPPFLAGS) -D PRINT_DEBUG=3
bin/debug: bin/headers/kernel/main.cl.h $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

bin/release: CCFLAGS  := $(CCFLAGS)  -O3
bin/release: LDFLAGS  := $(LDFLAGS)  -O3
bin/release: CPPFLAGS := $(CPPFLAGS) -D PRINT_DEBUG=0
bin/release: bin/headers/kernel/main.cl.h $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

bin/headers/%.cl.h: %.cl $(COMMON_OBJS)
	@mkdir -p $(dir $@)
	
	echo "char $(subst .,_,$(subst /,_,$<))_src[] = " > $@
	cat $< $(COMMON_SRCS) | $(CPP) $(CPPFLAGS) - | sed 's/\\/\\\\/g' | sed 's/\"/\\\"/g' | sed -e 's/\(.*\)/"\1\\n"/' >> $@
	echo ";" >> $@

bin/deps/%.cl.d: %.cl $(COMMON_OBJS)
	@mkdir -p $(dir $@)
	$(CPP) $(CPPFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%.d,bin/headers/%.h,$@)

bin/objs/%.c.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -o $@ $<

bin/deps/%.c.d: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CCFLAGS) -M -o $@ $< -MT $(patsubst bin/deps/%.d,bin/objs/%.o,$@)

