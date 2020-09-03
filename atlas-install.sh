./00boot
./configure --prefix=${VIRTUAL_ENV} --enable-swig-python --enable-lalinspiral --enable-lalburst --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps --enable-lalframe --enable-lalmetaio --disable-lalpulsar
make -j 16
make install
