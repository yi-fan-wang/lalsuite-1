installpath=/Users/wyf/Gwaves/lscsoft/opt/ecclal

./00boot
./configure --prefix=${installpath} --enable-swig-python --enable-lal --enable-lalsimulation --disable-lalburst --disable-lalinspiral --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps --enable-lalframe --disable-lalmetaio --disable-lalpulsar CFLAGS='-Wno-unused-but-set-variable -Wno-uninitialized -Wno-unused-parameter -Wno-sign-compare -Wno-unused-variable'
make -j 4
make install
