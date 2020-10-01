Notes on the usage of the Matlab interface to CUTEr
---------------------------------------------------

1. g77

To use MEX files with g77, make sure to load libg2c.so by setting PACKLIBS to
    
    -Lpath/to/g77/libs -lg2c

in $MYCUTER/bin/mx. Typically, libg2c.so is a system library located in a
directory such as /usr/lib so that there -lg2c should be sufficient.


2. gfortran

To use MEX files with gfortran, make sure to load libgfortran.so by setting
PACKLIBS to
    
    -Lpath/to/g77/libs -lgfortran

in $MYCUTER/bin/mx. Typically, libgfortran.so is a system library located in a
directory such as /usr/lib so that there -lgfortran should be sufficient.


3. g95

To use MEX files with g95, make sure to load libf95.a by setting PACKLIBS to

    -Lpath/to/g95/libs -lf95

in $MYCUTER/bin/mx, where path/to/g95/libs is the directory that contains
the libf95.a library shipped with g95. For instance, if you have unpacked
to g95 tar file in /home/myself/g95, the relevant directory is something
like

    /home/myself/g95/g95-install/lib/gcc-lib/i686-pc-linux-gnu/4.0.3
