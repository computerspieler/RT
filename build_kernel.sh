#!/bin/bash

make &&
gdb --args ./bin/debug kernel common
