#!/bin/bash

make debug &&
#gdb --args ./bin/debug scene/simple_plane.obj
gdb --args ./bin/debug scene/plane.obj
#scene/Internat.obj
#scene/bidirectional_rt.obj
#scene/Ball.obj
