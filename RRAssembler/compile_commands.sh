#!/bin/bash

cython3 best_path.pyx
g++ -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python3.10 -I/usr/local/lib/python3.10/dist-packages/numpy/core/include -o best_path.so best_path.c

cython3 dysgu_repeats.pyx
g++ -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python3.10 -o dysgu_repeats.so dysgu_repeats.c

cython3 recursive_find.pyx
g++ -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python3.10 -I/usr/local/lib/python3.10/dist-packages/pysam -I/usr/local/lib/python3.10/dist-packages/numpy/core/include -o recursive_find.so recursive_find.c

