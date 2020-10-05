#!/bin/bash

# Build script for GALAHAD 
# version for Bourne/bash shell

# syntax: install_galahad.sh

# N. Gould, D. Orban & Ph. Toint
# ( Last modified on 8 June 2008 at 07:00 GMT )

#  check input arguments (if any)

if [ $# != 0 ]; then
   echo "Use: install_galahad"
   exit 1
fi

#  function to create missing symblic links

galahad_create_missing_link () { 
 if [[ -f $1 && ! -L $2 ]] ;
   then echo "creating missing link $2" ;
#  else echo "link $2 already exists" ;
 fi ;
}

export ARCH=$PWD/arch/sh

#  determine the platform and operating system used

CORRECT_PLATFORM="false"
while [ $CORRECT_PLATFORM == "false" ]; do 
   echo ' Select platform'
   echo
   echo '   (1) Compaq (DEC) alpha'
   echo '   (2) Cray'
   echo '   (3) HP workstation'
   echo '   (4) IBM RS/6000'
   echo '   (5) PC'
   echo '   (6) PC with 64-bit Itanium processor'
   echo '   (7) PC with 64-bit Opteron processor'
   echo '   (8) PC with 64-bit Athlon processor'
   echo '   (9) SGI workstation'
   echo '  (10) SUN workstation'
   echo '  (11) MAC OS/X'   

   read CHOICE
   
   case  $CHOICE  in
       "1")
            CORRECT_PLATFORM="true"
            export MCH="alp"
            export MACHINE="Compaq Alpha"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) Digital Unix'
               echo '   (2) Tru-64 Unix'
               echo '   (3) Linux'
               
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="dux"
                  export OPSYS="Digital Unix"
;;
               2)
                  CORRECT_OS="true"
                  export OS="t64"
                  export OPSYS="Tru-64 Unix"
;;
               3)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 3'
               esac
            done
;;
       "2")
            CORRECT_PLATFORM="true"
            export MCH="cry"
            export MACHINE="CRAY T3E"
            export OS="unc"
            export OPSYS="UNICOS"
;;
       "3")
            CORRECT_PLATFORM="true"
            export MCH="hp"
            export MACHINE="HP workstation"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) HP-UX'
               echo '   (2) Linux'
            
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="hpu"
                  export OPSYS="HP-UX"
;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
;;
       "4")
            CORRECT_PLATFORM="true"
            export MCH="rs6"
            export MACHINE="IBM RS6000"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) AIX'
               echo '   (2) Linux'
           
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="aix"
                  export OPSYS="AIX"
;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
;;
       "5")
            CORRECT_PLATFORM="true"
            export MCH="pc"
            export MACHINE="Intel-like PC"

        CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) Windows 2000/XP with MinGW/Msys'
               echo '   (2) Linux'
           
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="mgw"
                  export OPSYS="MGW"
;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
;;
       "6")
            CORRECT_PLATFORM="true"
            export MCH="i64"
            export MACHINE="Intel-like PC with 64-bit Itanium processor"
            export OS="lnx"
            export OPSYS="Linux"
;;
       "7")
            CORRECT_PLATFORM="true"
            export MCH="opt"
            export MACHINE="AMD-like PC with 64-bit Opteron processor"
            export OS="lnx"
            export OPSYS="Linux"
;;
       "8")
            CORRECT_PLATFORM="true"
            export MCH="ath"
            export MACHINE="AMD-like PC with 64-bit Athlon processor"
            export OS="lnx"
            export OPSYS="Linux"
;;
       "9")
            CORRECT_PLATFORM="true"
            export MCH="sgi"
            export MACHINE="SGI workstation"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) IRIX'
               echo '   (2) Linux'
            
               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="irx"
                  export OPSYS="IRIX"
;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
;;
       "10")
            CORRECT_PLATFORM="true"
            export MCH="sun"
            export MACHINE="Sun workstation"

            CORRECT_OS="false"
            while [ $CORRECT_OS == "false" ]; do 

               echo ' Select operating system'
               echo
               echo '   (1) Solaris'
               echo '   (2) Linux'

               read CHOICE

               case  $CHOICE  in
               1)
                  CORRECT_OS="true"
                  export OS="sol"
                  export OPSYS="Solaris"
;;
               2)
                  CORRECT_OS="true"
                  export OS="lnx"
                  export OPSYS="Linux"
;;
               *)
                  echo ' Please give an integer between 1 and 2'
               esac
            done
;;
       "11")
             CORRECT_PLATFORM="true"
             export MCH="mac"
             export MACHINE="Mac"
             export OS="osx"
             export OPSYS="MacOSX"
;;
       *)
         echo ' Please give an integer between 1 and 11'
   esac
done

unalias source 2>/dev/null

source $ARCH/system.$OS

echo ' Select compiler'
echo

COMP=( `$LS $ARCH/compiler.${MCH}.${OS}.* $ARCH/compiler.all.all.*` )
NUMBER=${#COMP[@]}
LIST=( ${COMP[@]} )
let count=-1
for i  in  ${COMP[@]}; do
  (( count++ ))
  COMP[$count]="`$SED q $i | $SED 's/^[# ]*//'`"
done

CORRECT_COMPILER="false"

let count=-1
for i  in  ${LIST[@]}; do
  (( count++ ))
  let counter=count+1
  echo "        ($counter) ${COMP[$count]}"
done

while [ $CORRECT_COMPILER == "false" ]; do 
   read CHOICE
   let CHOICE=CHOICE-1
   if (( 0 <= CHOICE && CHOICE < NUMBER )); then
     CORRECT_COMPILER="true"
     COMPILER=${LIST[$CHOICE]##*/}
     CMP=${COMPILER##*\.} #${LIST[$CHOICE]%*\.}
     COMPUSED="${COMP[$CHOICE]}"
   else
     echo " Please give an integer between 1 and $NUMBER"
   fi
done

source $ARCH/$COMPILER

export GALAHAD=`dirs -l`
export GALAHAD=`echo $GALAHAD | $SED 's"/tmp_mnt""'`

VERSION=${MCH}.${OS}.${CMP}
#PREFIX=$VERSION
#OD='$GALAHAD'/objects/$VERSION

NEED_CUTER="false"
NEED_AMPL="false"
CORRECT_SUBSET="false"

if [[ $VERSION != "pc.lnx.frt" ]]; then
   while [ $CORRECT_SUBSET == "false" ]; do 
      echo ' Select subset of GALAHAD packages to be installed'
      echo ' (the chosen subset will optionally be installed below)'
      echo
      echo '     (1) Everything'
      echo '     (2) Everything for SIF/CUTEr'
      echo '     (3) Everything for AMPL'
      echo '     (4) LANCELOT B and its interface to SIF'
      echo '     (5) LANCELOT B and its interface to AMPL'
      echo '     (6) Just LANCELOT B'
      echo '     (7) The QP packages and their interfaces to CUTEr'
      echo '     (8) The QP packages and their interfaces to AMPL'
      echo '     (9) Just the QP packages and their dependencies'
      echo '    (10) FILTRANE and its interface to CUTEr'
      echo '    (11) FILTRANE and its interface to AMPL'
      echo '    (12) Just FILTRANE and its dependencies'
      
      read CHOICE
      
      case  $CHOICE  in
          "1")
               CORRECT_SUBSET="true"
               SUBSET="all"
               NEED_CUTER="true"
               NEED_AMPL="true"
;;
          "2")
               CORRECT_SUBSET="true"
               SUBSET="all_cuter"
               NEED_CUTER="true"
;;
          "3")
               CORRECT_SUBSET="true"
               SUBSET="all_ampl"
               NEED_AMPL="true"
;;
          "4")
               CORRECT_SUBSET="true"
               SUBSET="lanb_sif"
;;
          "5")
               CORRECT_SUBSET="true"
               SUBSET="lanb_ampl"
               NEED_AMPL="true"
;;
          "6")
               CORRECT_SUBSET="true"
               SUBSET="lanb"
;;
          "7")
               CORRECT_SUBSET="true"
               SUBSET="qps_cuter"
               NEED_CUTER="true"
;;
          "8")
               CORRECT_SUBSET="true"
               SUBSET="qps_ampl"
               NEED_AMPL="true"
;;
          "9")
               CORRECT_SUBSET="true"
               SUBSET="qps"
;;
          "10")
               CORRECT_SUBSET="true"
               SUBSET="filtrane_cuter"
               NEED_CUTER="true"
;;
          "11")
               CORRECT_SUBSET="true"
               SUBSET="filtrane_ampl"
               NEED_AMPL="true"
;;
          "12")
               CORRECT_SUBSET="true"
               SUBSET="filtrane"
;;
          *)
            echo ' Please give an integer between 1 and 9'
      esac
   done
else
   while [[ $CORRECT_SUBSET == "false" ]]; do 
      echo ' Select subset of GALAHAD packages to be installed'
      echo ' (the chosen subset will optionally be installed below)'
      echo
      echo '    (1) Everything for SIF/CUTEr'
      echo '    (2) LANCELOT B and its interface to SIF'
      echo '    (3) Just LANCELOT B'
      echo '    (4) The QP packages and their interfaces to CUTEr'
      echo '    (5) Just the QP packages and their dependencies'
      echo '    (6) FILTRANE and its interface to CUTEr'
      echo '    (7) Just FILTRANE and its dependencies'
      
      read CHOICE
      
      case  $CHOICE  in
          "1")
               CORRECT_SUBSET="true"
               SUBSET="all_cuter"
               NEED_CUTER="true"
;;
          "2")
               CORRECT_SUBSET="true"
               SUBSET="lanb_sif"
;;
          "3")
               CORRECT_SUBSET="true"
               SUBSET="lanb"
;;
          "4")
               CORRECT_SUBSET="true"
               SUBSET="qps_cuter"
               NEED_CUTER="true"
;;
          "5")
               CORRECT_SUBSET="true"
               SUBSET="qps"
;;
          "6")
               CORRECT_SUBSET="true"
               SUBSET="filtrane_cuter"
               NEED_CUTER="true"
;;
          "7")
               CORRECT_SUBSET="true"
               SUBSET="filtrane"
;;
          *)
            echo ' Please give an integer between 1 and 5'
      esac
   done
fi

#  check to see which CUTEr is required

if [[ $NEED_CUTER == "true" ]]; then
    CUTERUSED=$MYCUTER
    
    YESNO=""
    
    echo ' By default, the CUTEr you wish to use is installed in'
    echo "  $CUTERUSED"
    echo ' Is this OK (Y/n)?'
    while [[ $YESNO != 'Y' && $YESNO != 'N' ]]; do
        read YESNO
        if [[ $YESNO == "" ]]; then
        YESNO="Y"
    fi
        YESNO=`echo $YESNO | tr a-z A-Z`
    done
    
    CUTER_FOUND='no'

    while [[ $CUTER_FOUND == 'no' ]]; do
       if [[ $YESNO == 'N' ]]; then
           echo " Enter alternative directory for CUTEr:"
           read CUTERUSED
       fi
       if [[ -e $CUTERUSED ]]; then
          CUTER_FOUND='yes'
       else
          echo " The CUTEr directory you have specified does not exist."
          YESNO='N'
       fi
    done
else
    CUTERUSED=
fi

#  check to see if AMPLDIR is defined

if [[ $NEED_AMPL == "true" ]]; then
    if [[ ${AMPLDIR+set} == 'set' ]]; then
      AMPLLIBDIR=$AMPLDIR
      YESNO=""
      echo ' By default, the AMPL interface library you wish to use is in'
      echo "  $AMPLLIBDIR"
      echo ' Is this OK (Y/n)?'
      while [[ $YESNO != 'Y' && $YESNO != 'N' ]]; do
          read YESNO
          if [[ $YESNO == "" ]]; then
        YESNO="Y"
      fi
          YESNO=`echo $YESNO | tr a-z A-Z`
      done
    else
       echo ' You plan to use the AMPL interface but the'
       echo ' AMPLDIR environment variable is not currently set'
       YESNO='N'
    fi
    
    CORRECT_AMPLLIBDIR="false"
    while [[ $CORRECT_AMPLLIBDIR == "false" ]]; do 
       if [[ $YESNO == 'N' ]]; then
          echo ' Please give the name of the directory'
          echo ' containing the AMPL interface library:'
          read AMPLLIBDIR
       fi
       if [[ -e $AMPLLIBDIR/amplsolver.a ]]; then
          CORRECT_AMPLLIBDIR="true"
       else
          echo " The directory $AMPLLIBDIR"
          echo " does not appear to contain a working AMPL interface library."
          YESNO='N'
       fi
    done

    echo ' Select C compiler'
    echo

    CCOMP=( `$LS $ARCH/ccompiler.${MCH}.${OS}.* $ARCH/ccompiler.all.all.*`)
    NUMBER=${#CCOMP[@]}
    LIST=( ${CCOMP[@]} )
    let count=-1
    for i  in  ${CCOMP[@]}; do
      (( count++ ))
      CCOMP[$count]="`$SED q $i | $SED 's/^[# ]*//'`"
    done

    CORRECT_CCOMPILER="false"
    while [[ $CORRECT_CCOMPILER == "false" ]]; do 

        let count=-1
        for i  in  ${LIST[@]}; do
          (( count++ ))
          echo "        ($count) ${CCOMP[$count]}"
        done

        read CHOICE
        
        i=0
        while [[ $i < $NUMBER &&  $CORRECT_CCOMPILER == "false" ]]; do
           (( i++ ))
           if [[ $CHOICE == $i ]]; then
             CORRECT_CCOMPILER="true"
             CHOICE=$i
           fi
        done
        if [[ $CORRECT_CCOMPILER == "true" ]]; then
          CCOMPILER=${LIST[$CHOICE]##*/}
          CMP=${CCOMPILER##*\.} #${LIST[$CHOICE]%*\.}
          #CCOMPILER=$LIST[$CHOICE]:t
          #CMP=$LIST[$CHOICE]:e
          CCOMPUSED="${CCOMP[$CHOICE]}"
        else
          echo " Please give an integer between 1 and $NUMBER"
        fi
    done
else
    AMPLLIBDIR=
    CCOMPILER=ccompiler.all.all.gcc
fi

source $ARCH/$CCOMPILER

#  create architecture-dependent object and module directories

OBJDIR=$GALAHAD/objects/$VERSION
MODDIR=$GALAHAD/modules/$VERSION

echo "$MACHINE ($OPSYS) $COMPUSED" > $GALAHAD/versions/$VERSION

if [[ ! -e $OBJDIR ]]; then
    $MKDIR $OBJDIR 
    $MKDIR $OBJDIR/double $OBJDIR/single
else
    if [[ ! -e $OBJDIR/double ]]; then
    $MKDIR $OBJDIR/double
    fi
    if [[ ! -e $OBJDIR/single ]]; then
    $MKDIR $OBJDIR/single
    fi
fi

if [[ ! -e $MODDIR ]]; then
    $MKDIR $MODDIR 
    $MKDIR $MODDIR/double $MODDIR/single
else
    if [[ ! -e $MODDIR/double ]]; then
    $MKDIR $MODDIR/double
    fi
    if [[ ! -e $MODDIR/single ]]; then
    $MKDIR $MODDIR/single
    fi
fi

#  for compilers based on the Edinburgh portable compiler, set
#  the work.pcl file

if [[ $CMP == "e90" ]]; then

    echo "work.pc" > $MODDIR/double/work.pcl
    echo "$OBJDIR/double/work.pc" >> $MODDIR/double/work.pcl
    epcfcem<<!
    cr $OBJDIR/double/work.pc
!

    echo "work.pc" > $MODDIR/single/work.pcl
    echo "$OBJDIR/single/work.pc" >> $MODDIR/single/work.pcl
    epcfcem<<!
    cr $OBJDIR/single/work.pc
!

#elif [[ $CMP == "ifc" ]]; then

#    echo "work.pc" > $MODDIR/double/work.pcl
#    echo "$OBJDIR/double/work.pc" >> $MODDIR/double/work.pcl
#    ifccem<<!
#    cr $OBJDIR/double/work.pc
#!

#    echo "work.pc" > $MODDIR/single/work.pcl
#    echo "$OBJDIR/single/work.pc" >> $MODDIR/single/work.pcl
#    ifccem<<!
#    cr $OBJDIR/single/work.pc
#!

fi

#  write out the galahad/bin/sys file for this architecture

SYSFILE=$GALAHAD/bin/sys/$VERSION

echo 'RM="'$RM'"'                                                  >  $SYSFILE
echo 'MAKE="'$MAKE'"'                                              >> $SYSFILE
echo 'CAT="'$CAT'"'                                                >> $SYSFILE
echo 'SED="'$SED'"'                                                >> $SYSFILE
echo 'MV="'$MV'"'                                                  >> $SYSFILE
echo 'LS="'$LS'"'                                                  >> $SYSFILE
echo 'FORTRAN="'$FORTRAN'"'                                        >> $SYSFILE
MOD='$GALAHAD/modules/'$VERSION'/$PRECIS'
FFLAGS="$LIBCMD"' '`eval echo $MODCMD`' '"$F90"
echo 'FFLAGS="'$FFLAGS'"'                                          >> $SYSFILE
echo 'PROBFLAGS="'$FFLAGS' '$BASIC' '$OPTIMIZATION' '$F77'"'       >> $SYSFILE
echo 'CUTERUSED="'$CUTERUSED'"'                                    >> $SYSFILE
echo 'BLAS="'$BLAS'"'                                              >> $SYSFILE
echo 'LAPACK="'$LAPACK'"'                                          >> $SYSFILE
echo 'HSL="'$HSL'"'                                                >> $SYSFILE
echo 'METIS="'$METIS'"'                                            >> $SYSFILE

#  write out the galahad/makefile/ file for this architecture

MAKEFILE=$GALAHAD/makefiles/$VERSION

echo ' '                                                           >  $MAKEFILE
echo '#  Architecture dependent makefile'                          >> $MAKEFILE
echo '#  (automatically generated by install_galahad)'             >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'VERSION = '$VERSION                                          >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Basic system commands'                                    >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'CP = '$CP                                                    >> $MAKEFILE
echo 'MV = '$MV                                                    >> $MAKEFILE
echo 'RM = '$RM                                                    >> $MAKEFILE
echo 'SED = '$SED                                                  >> $MAKEFILE
echo 'GREP = '$GREP                                                >> $MAKEFILE
echo 'AR = '$AR                                                    >> $MAKEFILE
echo 'RANLIB = '$RANLIB                                            >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Directory for binaries'                                   >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'PRECIS = double'                                             >> $MAKEFILE
echo 'OBJ = $(GALAHAD)/objects/$(VERSION)/$(PRECIS)'               >> $MAKEFILE
echo 'OBJS = $(GALAHAD)/objects/$(VERSION)/single'                 >> $MAKEFILE
echo 'OBJD = $(GALAHAD)/objects/$(VERSION)/double'                 >> $MAKEFILE
echo 'MOD = $(GALAHAD)/modules/$(VERSION)/$(PRECIS)'               >> $MAKEFILE
echo 'SEDS = $(GALAHAD)/seds/$(PRECIS).sed'                        >> $MAKEFILE
echo 'MVMODS = '"$MVMODS"                                          >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Compiler options'                                         >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'FORTRAN = '$FORTRAN                                          >> $MAKEFILE
echo 'BASIC = '$BASIC                                              >> $MAKEFILE
MODTMP="$LIBCMD"' '`echo $MODCMD | $SED 's/MOD/(MOD)/g'`
echo 'MODULES = '$MODTMP                                           >> $MAKEFILE
echo 'OPTIMIZATION = '$OPTIMIZATION                                >> $MAKEFILE
echo 'NOOPTIMIZATION = '$NOOPTIMIZATION                            >> $MAKEFILE
echo 'DEBUG = '$DEBUG                                              >> $MAKEFILE
echo 'F77 = '$F77                                                  >> $MAKEFILE
echo 'F90 = '$F90                                                  >> $MAKEFILE
echo 'F95 = '$F95                                                  >> $MAKEFILE
echo 'NOFMAIN = '$NOFMAIN                                          >> $MAKEFILE
echo 'USUAL = '$USUAL                                              >> $MAKEFILE
echo 'SPECIAL = '$SPECIAL                                          >> $MAKEFILE
echo 'F77SUFFIX = '$F77SUFFIX                                      >> $MAKEFILE
echo 'F95SUFFIX  = '$F95SUFFIX                                     >> $MAKEFILE
echo 'TIMER = '$TIMER                                              >> $MAKEFILE
echo 'NOT95 = '$NOT95                                              >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'AMPLDIR   = '$AMPLLIBDIR                                     >> $MAKEFILE
echo 'CC        = '$CC                                             >> $MAKEFILE
echo 'CCBASIC   = '$CCBASIC                                        >> $MAKEFILE
echo 'CCONDEF   = '$CCONDEF                                        >> $MAKEFILE
echo 'CCDEBUG   = '$CCDEBUG                                        >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Libraries'                                                >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'BLAS = '$BLAS                                                >> $MAKEFILE
echo 'LAPACK = '$LAPACK                                            >> $MAKEFILE
echo 'HSL = '$HSL                                                  >> $MAKEFILE
echo 'METIS = '$METIS                                              >> $MAKEFILE
echo 'CUTERUSED = '$CUTERUSED                                      >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Shell used'                                               >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'BINSHELL = '$BINSHELL                                        >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Set directories for optional packages'                    >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'include $(GALAHAD)/src/makedefs/packages'                    >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo '#  Body of makefile'                                         >> $MAKEFILE
echo ' '                                                           >> $MAKEFILE
echo 'include $(PWD)/makemaster'                                   >> $MAKEFILE

#  check that required symbolic links are in place

galahad_create_missing_link $GALAHAD/python/README \
                            $GALAHAD/doc/README.gui
galahad_create_missing_link $GALAHAD/doc/README \
                            $GALAHAD/README
galahad_create_missing_link $GALAHAD/python/galahad \
                            $GALAHAD/bin/galahad
galahad_create_missing_link $GALAHAD/bin/install_galahad \
                            $GALAHAD/install_galahad
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/dum/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/spec/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/string/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/lapack/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/sort/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/rand/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/scu/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/sils/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/smt/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/sym/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/auxiliary/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/ptrans/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/space/makemaster 
galahad_create_missing_link $GALAHAD/src/general/makemaster \
                            $GALAHAD/src/trans/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/gltr/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/lsqp/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/pre/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpa/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpb/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpc/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpd/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpp/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qpt/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/qtrans/makemaster 
galahad_create_missing_link $GALAHAD/src/qp/makemaster \
                            $GALAHAD/src/roots/makemaster 
galahad_create_missing_link $GALAHAD/src/lpsqp/makemaster \
                            $GALAHAD/src/lpqp/makemaster 
galahad_create_missing_link $GALAHAD/src/lpsqp/makemaster \
                            $GALAHAD/src/lpqpa/makemaster 
galahad_create_missing_link $GALAHAD/src/lpsqp/makemaster \
                            $GALAHAD/src/lpqpb/makemaster 
galahad_create_missing_link $GALAHAD/src/filtrane/makemaster \
                            $GALAHAD/src/nlpt/makemaster 
galahad_create_missing_link $GALAHAD/src/fastr/makemaster \
                            $GALAHAD/src/pqp/makemaster 
galahad_create_missing_link $GALAHAD/src/eqp/makemaster \
                            $GALAHAD/src/sbls/makemaster 
galahad_create_missing_link $GALAHAD/doc/README.trouble \
                            $GALAHAD/README.trouble
galahad_create_missing_link $GALAHAD/bin/uninstall_galahad \
                            $GALAHAD/uninstall_galahad
galahad_create_missing_link $GALAHAD/doc/README.SIF \
                            $GALAHAD/README.SIF
galahad_create_missing_link $GALAHAD/doc/README.windows \
                            $GALAHAD/README.windows

#  optionally compile the selected packages

YESNO=""
echo ' '
echo ' Do you now wish to compile the package subset you selected earlier (Y/n)?'
while [[ $YESNO != 'Y' && $YESNO != 'N' ]]; do
    read YESNO
    if [[ $YESNO == "" ]]; then
    YESNO="Y"
    fi
    YESNO=`echo $YESNO | tr a-z A-Z`
done

if [[ $YESNO == 'Y' ]]; then

    PREC=""
    echo ' '
    echo ' The package subset may be installed in either single or double precision'
    echo ' Which precision do you require for the installed subset ? '
    echo ' D for double precision (the default), S for single precision'
    while [[ $PREC != 'S' && $PREC != 'D' ]]; do
        read PREC
        if [[ $PREC == "" ]]; then
        PREC="D"
    fi
        PREC=`echo $PREC | tr a-z A-Z`
    done 
    if [[ $PREC == 'S' ]]; then
      PREC="single"
    else
      PREC="double"
    fi

    cd $GALAHAD/src/

    echo ' '
    echo "Installing the $PREC precision version"
    OPTIONS="-s -f $GALAHAD/makefiles/$VERSION"
    MACROS="PRECIS=$PREC PWD=$GALAHAD/src GALAHAD=$GALAHAD"
    case  $SUBSET  in
        "all")
            $MAKE $OPTIONS all $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "all_cuter")
            $MAKE $OPTIONS all_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "all_ampl")
            $MAKE $OPTIONS all_ampl $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "lanb_sif")
            $MAKE $OPTIONS lancelotb_sif $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "lanb_ampl")
            $MAKE $OPTIONS lancelotb_ampl $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "lanb")
            $MAKE $OPTIONS lancelotb $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "qps_cuter")
            $MAKE $OPTIONS qp_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "qps_ampl")
            $MAKE $OPTIONS qp_ampl $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "qps")
            $MAKE $OPTIONS qp $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "filtrane_cuter")
            $MAKE $OPTIONS filtrane_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "filtrane_ampl")
            $MAKE $OPTIONS filtrane_ampl $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "filtrane")
            $MAKE $OPTIONS filtrane $MACROS
            [[ $? != 0 ]] && exit 2
;;
    esac

#  optionally compile the selected packages in the other precision

    if [[ $PREC == 'single' ]]; then
      PREC="double"
    else
      PREC="single"
    fi

    YESNO2=""
    echo ' '
    echo "Do you also wish to install the $PREC precision version ? (N/y)"

    while [[ $YESNO2 != 'Y' && $YESNO2 != 'N' ]]; do
        read YESNO2
        if [[ $YESNO2 == "" ]]; then
        YESNO2="N"
    fi
        YESNO2=`echo $YESNO2 | tr a-z A-Z`
    done
    
    if [[ $YESNO2 == 'Y' ]]; then

        echo ' '
        echo "Installing the $PREC precision version"
        MACROS="PRECIS=$PREC PWD=$GALAHAD/src GALAHAD=$GALAHAD"
        case  $SUBSET  in
        "all")
            $MAKE $OPTIONS all_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "all_cuter")
            $MAKE $OPTIONS all_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "lanb_sif")
            $MAKE $OPTIONS lancelotb_sif $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "lanb")
            $MAKE $OPTIONS lancelotb $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "qps_cuter")
            $MAKE $OPTIONS qp_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "qps")
            $MAKE $OPTIONS qp $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "filtrane_cuter")
            $MAKE $OPTIONS filtrane_cuter $MACROS
            [[ $? != 0 ]] && exit 2
;;
        "filtrane")
            $MAKE $OPTIONS filtrane $MACROS
            [[ $? != 0 ]] && exit 2
;;
        esac

    fi
fi

echo ''
echo " Remember to set the environment variable"
echo "  GALAHAD to $GALAHAD"
echo " In addition, please update your MANPATH to include"
echo "    $GALAHAD/man"
echo " and your PATH to include"
echo "    $GALAHAD/bin"
echo ''
echo " (see $GALAHAD/doc/README.cshrc"
echo "  and $GALAHAD/doc/README.bashrc"
echo "  for examples on how to do this)"
echo ''

exit 0

#  function to create missing symblic links

galahad_create_missing_link () { 
 if [[ -f $1 && ! -L $2 ]] ;
   then echo "creating missing link $2" ;
#  else echo "link $2 already exists" ;
 fi ;
}

# N. Gould, D. Orban and Ph.L. Toint, 17th March, 2002.
# This version 31st October, 2005.

