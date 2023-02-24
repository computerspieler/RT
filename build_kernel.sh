#!/bin/bash

make debug &&
gdb --args ./bin/debug scene/plane.obj
#scene/Internat.obj
#scene/bidirectional_rt.obj
