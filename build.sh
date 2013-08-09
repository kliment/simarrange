#!/bin/sh
gcc -O3 -I . -Wno-unused-result -fno-builtin-strlen admesh/*.c simarrange.c -o ./sa -lm -lcv -lcxcore -lhighgui -largtable2
