# Springs -- A simple integrator for the motion of coupled and dampened springs #

Compilation: set CXX in the Makefile to your preferred C++ compiler, then compile with `make`.

Usage: `./springs -m M -k K -b B`, where `M` is the (shared) node mass, `K` the spring constant, and `B` the spring dampening coefficient. After running, the code will automatically invoke Python for output data visualization.
