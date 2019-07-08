installpath=/home/astro/Gwaves/lscsoft/opt/mylalfork-cor-tgr

./00boot
./configure --prefix=${installpath} --enable-swig-python --disable-lalstochastic --disable-lalxml --disable-lalinference --disable-laldetchar --disable-lalapps --disable-lalframe --disable-lalmetaio --disable-lalpulsar
make -j 8
make install
