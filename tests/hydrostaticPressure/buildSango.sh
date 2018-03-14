








#!/bin/bash
CC=icc CFLAGS="-O2 -mtune=native -mkl=sequential" POPT_CFLAGS="/usr/include" POPT_LIBS="/usr/lib64/libpopt.so" ./configure
make -j

