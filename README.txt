% NL-Bayes video denoising.

# ABOUT

* Author    : Pablo Arias <pariasm@gmail.com>
* Copyright : (C) 2015 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the NL-Bayes video denoising.

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS (not tested!). 

- Compilation. 
Automated compilation requires the make program.

- Library. 
This code requires the libpng library.

- Image format. 
Only the PNG format is supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Configure and compile the source code using cmake and make. 
It is recommended that you create a folder for building:

$ mkdir build; cd build
$ cmake ..
$ make

Binaries will be created in build/bin folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). It can 
be disabled by editing the root CMakeLists.txt file.

3. Run NL-Bayes video denoising. 
Type bin/vnlbayes --help for usage instructions.

4. Results are available in the file "measures.txt".

# ABOUT THIS FILE

Copyright 2015 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
