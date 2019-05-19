#ifndef GT_CONFIG_H
#define GT_CONFIG_H
#define GT_CC "cc (Ubuntu 4.9.1-16ubuntu6) 4.9.1"
#define GT_CFLAGS " -g -Wall -Wunused-parameter -pipe -fPIC -Wpointer-arith -Wno-unknown-pragmas -O3 -m32 -m64 -Werror"
#define GT_CPPFLAGS "-fno-stack-protector -U_FORTIFY_SOURCE -D_GNU_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHAVE_MEMMOVE -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -DLUA_DL_DLOPEN -DLUA_USE_MKSTEMP -DWITHOUT_CAIRO -DSQLITE_THREADSAFE=0 -DHAVE_SQLITE -I/vagrant/src -I/vagrant/obj -I/vagrant/src/external/zlib-1.2.8 -I/vagrant/src/external/md5-1.2/src -I/vagrant/src/external/lua-5.1.5/src -I/vagrant/src/external/luafilesystem-1.5.0/src -I/vagrant/src/external/lpeg-0.10.2 -I/vagrant/src/external/expat-2.0.1/lib -I/vagrant/src/external/bzip2-1.0.6 -I/vagrant/src/external/samtools-0.1.18 -I/vagrant/src/external/sqlite-3.8.7.1 -I/vagrant/src/external/tre/include -I/vagrant/src/external/sqlite-3.8.7.1"
#define GT_VERSION "1.5.10"
#define GT_MAJOR_VERSION 1
#define GT_MINOR_VERSION 5
#define GT_MICRO_VERSION 10
#endif
