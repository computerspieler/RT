#!/bin/bash

make &&
gdb --args ./bin/debug kernel/main.c common
