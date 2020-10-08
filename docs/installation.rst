.. _installation:

******************
Installation guide
******************

.. _sys-requirements:

System requirements
===================

Computers running the DYNAMITE models must meet the following hardware requirements.

**Minimum Hardware Requirements**

* Processor:
* Processor speed:
* Random access memory (RAM):
* GPU
* Hard disk capacity:

**Recommended Hardware Requirements**

* Processor:
* Processor speed:
* access memory (RAM):
* GPU
* Hard disk capacity:


.. _software-requirements:

Software requirements
=====================

The current version of DYNAMITE is based on routines which were written in Fortran. It is therefore necessary to have a pre-installed Fortran compiler. The best-suited Fortran compiler differs with different operation systems.

If you are using Linux, you will need ``ifort`` or ``gcc-gfortran``.

If you are using macOS, you will need ``gfortran``, which can be installed in a number of different ways. In the following, we explain the installation via Homebrew and via MacPorts. Homebrew can be used to install the latest GCC and all additional libraries in the following way::

    brew update
    brew install gcc

We can check if gfortran is installed by typing in the Terminal::

    man -k fortran

with the output: ``gfortran(1) - GNU Fortran compiler``. The command ::

    which gfortran

returns its location, for example: ``/usr/local/bin/gfortran``. Alternatively, you can also install GCC with MacPorts. The installation is then a bit different::

    sudo port selfupdate
    sudo upgrade outdated
    port search --name --line --glob 'gcc*'

lists the available GCC version. For example, having as the latest one gcc9, one should then proceed as follows::

    sudo port install gcc9
    sudo port select --set gcc mp-gcc9

to install it and make it the default version (including gfortran). Check and display the path with::

    gfortran -v
    ...
    which gfortran

which returns something like ``/opt/local/bin/gfortran``.

In any case, when using macOS/Linux, the basic requirements for this code are reasonably current versions of python (tested with Legacy Python 2.7, Python 3.6, Python 3.7, and Python 3.8) and the following python libraries::

  numpy
  glob
  re
  os
  matplotlib
  astropy
  math
  shutil
  scipy
  ast
  argparse
  configparser

Pre-Installation Checklist
--------------------------

...

Known problems
--------------

Make sure that the compiler is installed with the proper libraries.


.. _install-procedure:

Installation and Configure Procedure
====================================

Download from `github <https://github.com/dynamics-of-stellar-systems/triaxschwarz>`_, unzip and move the DYNAMITE code to the directory in which you want to install it. Make sure that your system fulfills the software requirements specified in Section~\ref{soft_requirements} (in particular the Fortran compiler).


Installation of GALAHAD
-----------------------

GALAHAD is a "library of thread-safe Fortran 90 packages for large-scale nonlinear optimization". The DYNAMITE code comes with Version 2.3.  An updated version of GALAHAD could be obtained `here <http://www.galahad.rl.ac.uk/doc.html>`_ (last updated in 2018), but the most recent version seems to not work. The GALAHAD package included in DYNAMITE can be found in the folder ``.../legacy_fortran``.

For the installation go into the folder ``.../legacy_fortran/galahad-2.3/`` and type ::

    ./install_galahad

In the following installation, a number of prompts start. The answers differ for the different operation system and are shown in the following.


Install Galahad, Version 2.3 - Prompt answers for Ubuntu 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prompts from ``./install_galahad``. The answers for the recommended installation are marked in bold.

**Select platform**


* Compaq (DEC) alpha
* Cray
* HP Workstation
* IBM RS/6000
* **PC**
* ...


**Select operating system**

* Windows 2000/XP with MinGW/Msys
* **Linux**

**Select compiler**

When using Ubuntu:

* Windows 2000/XP with MinGW/Msys
* **Linux**


**Select subset of GALAHAD packages to be installed (the chosen subset will optionally be installed below)**

* Everything
* Everything for SIF/CUTEr
* Everything for AMPL
* LANCELOT B and its interface to SIF
* LANCELOT B and its interface to AMPL
* Just LANCELOT B
* **The QP packages and their interfaces to CUTEr**
* The QP packages and their interfaces to AMPL
* Just the QP packages and their dependencies
* FILTRANE and its interface to CUTEr
* FILTRANE and its interface to AMPL
* Just FILTRANE and its dependencies

**By default, the CUTEr you wish to use is installed in**

* y(es)
* **n(o)**

**Enter alternative directory for CUTEr:**

  | ``/Users/.../dynamite/legacy_fortran/cuter`` (Note: Put your full directory path here)

**Do you now wish to compile the package subset you selected earlier?**

* **y(es)**
* n(o)

**The package subset may be installed in either single or double precision. Which precision do you require for the installed subset?**

* **D for double precision**
* S for single precision

**Do you also wish to install the single precision version?**

* y(es)
* **n(o)**

Install Galahad, Version 2.3 - Prompt answers for MacOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Prompts from ``./install_galahad``. The answers for the recommended installation are marked in bold.

**Select platform**

* Compaq (DEC) alpha
* Cray
* HP Workstation
* IBM RS/6000
* PC
* PC with ..
* PC with
* PC with
* SGI workstation
* SUN workstation
* **MAC OS/X**

**Select compiler**

When using MacOS:

* NAG f90
* NAG f95
* AbSoft f95
* GNU g95 under OS/X
* **GNU gfortran under OS/X**
* Intel ifort (previously ifc) under Mac OsX

**Select subset of GALAHAD packages to be installed (the chosen subset will optionally be installed below)**

* Everything
* Everything for SIF/CUTEr
* Everything for AMPL
* LANCELOT B and its interface to SIF
* LANCELOT B and its interface to AMPL
* Just LANCELOT B
* **The QP packages and their interfaces to CUTEr**
* The QP packages and their interfaces to AMPL
* Just the QP packages and their dependencies
* FILTRANE and its interface to CUTEr
* FILTRANE and its interface to AMPL
* Just FILTRANE and its dependencies

**By default, the CUTEr you wish to use is installed in**

* y(es)
* **n(o)**

**Enter alternative directory for CUTEr:**

  | ``/Users/.../dynamite/legacy_fortran/cuter`` (Note: Put your full directory path here)

**Do you now wish to compile the package subset you selected earlier?**

* **y(es)**
* n(o)

**The package subset may be installed in either single or double precision. Which precision do you require for the installed subset?**

* **D for double precision**
* S for single precision

**Do you also wish to install the single precision version?**

* y(es)
* **n(o)**


Finalizing the installation of GALAHAD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set environment variables and path as prompted at the end of successful Galahad installation e.g. in your .bashrc file.

  | COMMENT: At the end of the installation, the output hints towards the README.bashrc and README.cshrc files, to look up how to set the environment variables correctly. The content of these files however is a bit confusing, so maybe this could be changed.

**Example: GALAHAD environment variables**

Output from GALAHAD::

    Remember to set the environment variable
    GALAHAD to /home/fs71474/sthater/triaxschwarz/galahad-2.3
    In addition, please update your MANPATH to include
    /home/fs71474/sthater/triaxschwarz/galahad-2.3/man
    and your PATH to include
    /home/fs71474/sthater/triaxschwarz/galahad-2.3/bin

Update in .bashrc::

    export GALAHAD="/home/fs71474/sthater/triaxschwarz/galahad-2.3/galahad"
    export MANPATH="$/home/fs71474/sthater/triaxschwarz/galahad-2.3/man:$/home/fs71474/sthater/triaxschwarz/galahad-2.3/man"
    export PATH="$/home/fs71474/sthater/triaxschwarz/galahad-2.3/bin:$PATH"


Compiling the DYNAMITE code
---------------------------

Go back to ``.../triaxschwarzschild``. Before you proceed, it is necessary to make three changes to the ``Makefile``:

* Change the local path of Galahad (``GALAHADDIR``) to something like ``/Users/.../triaxschwarz/triaxschwarzschild/galahad-2.3``.
* Select the appropriate choice of ``GALAHADTYPE`` variable depending on your system (possible options are commented out)
* Look for the definition of the ``all:`` (this should be right after the definition of the ``GALAHADTYPE`` variable). Make sure that ``triaxgasnnls`` is **NOT** in the list.

If you install and run DYNAMITE on your own computer, there seems to be a memory allocation problem when building the orbit library. This problem has currently been tackled by adding a line in the source code (the line ``print*, t1,t2,t3`` right after ``losvel(:) = t1 * vel(:,1) + t2 * vel(:,2) + t3 * vel(:,3)`` in the subroutine ``project_n(type,pos,vel,proj,losvel,n)``).

Proceed with the following command from the terminal::

    make all

DYNAMITE should now be installed and ready to be run! You can try the test run as explained in Test_ run.


Uninstalling DYNAMITE from the system
-------------------------------------
With the following command from the terminal::

    make distclean

all compiled Fortran codes are removed.



Post-Installation
=================

Post-installation checklist
---------------------------

* Check that all files are there
* Check which NNLS you want to be used

.. _Test:

Test run
--------

You can have a test run of the Dynamite code and the analysis scripts on the S0 galaxy NGC 6278. In the end you should get similar plots to the ones shown in `Zhu et al. 2018, MNRAS, 473, 3000 <https://arxiv.org/pdf/1709.06649.pdf>`_.
For this testrun, we have created a data directory in ``/Users/.../triaxschwarz/model_example/NGC6278``, containing all the necessary data. This directory only includes the folder ``infil``, which contains the input files of the DYNAMITE code. If you run the code with your own data, make sure that your galaxy folder (named by the object name) has all input files with the parameters set properly for your galaxy.
Before starting the DYNAMITE run, you need to change the directories in the ``my_config.ini``-file which is located in the folder ``/Users/.../triaxschwarz/schwpy``.

This directory specifies the folder of your Dynamite code (the Fortran scripts)::

    exepath = '/Users/.../triaxschwarz/triaxschwarzschild/'

This directory specifies the parent folder of your future models::
    w_dir="/Users/sabine/Triaxschwarzschild/model_example/"

``object="NGC6278"``
Here you declare the name of the galaxy that you want to model. The name should also be used for the data folder which contains your galaxy input files.

Finally, we can start the first models now. In the terminal, type ::

    python schw_run.py -f my_config.ini
    python schw_run.py -f -m my_config.ini

Note: When running the code, the NNLS skips one M/L and the model for M/L=0.8 is not made. Only when I delete all M/L files and start ``schw_run`` again, it makes all models. ::

    python schw_run.py -d my_config.ini
    python schw_run.py -d -m my_config.ini
    python schw_run.py -p  my_config.ini


Troubleshooting
===============

...
