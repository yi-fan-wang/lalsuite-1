lalpath=/home/astro/Gwaves/lscsoft/opt/master/include
lallibpath=/home/astro/Gwaves/lscsoft/opt/master/lib

gcc -pthread -I${lalpath} GRFlagsTest.c -L${lallibpath} -lm -llal -llalsimulation

