#!/usr/bin/env bash
# Builds the voter-model executable v2d (same flags as the original Makefile: -O3, -DDUMP
# so that every realization is written to an AConf_* file).
cd "$(dirname "$0")" && cc -O3 -DDUMP -o v2d voter2d.c evolution.c lookuptable.c moving.c statistics.c rnd.c debug.c -lm
