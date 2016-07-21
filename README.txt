NON LOCAL BAYESIAN IMAGE/VIDEO DENOISING.

=====
ABOUT
=====

* Author    : Pablo Arias <pariasm@gmail.com>
* Copyright : (C) 2015 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

========
OVERVIEW
========

This source code provides an implementations of some NL Bayesian image and 
video denoising algorithms.

=========================
UNIX/LINUX/MAC USER GUIDE
=========================

The code is compilable on Unix/Linux and hopefully on Mac OS (not tested!). 

Compilation: requires the cmake and make programs.

Dependencies: CBLAS, LAPACKE, OpenMP [can be disabled], iio.
 
-----
Usage
-----

1/
Download the code package and extract it. Go to that directory. 

2/
Configure and compile the source code using cmake and make. 
It is recommended that you create a folder for building:

$ mkdir build; cd build
$ cmake ..
$ make

Binaries will be created in build/bin folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). It can 
be disabled by editing the root CMakeLists.txt file.

3/
Run NL-Bayes video denoising. 
To denoise video

bin/vnlbayes \
-i path/to/sequence/frame_%03d.png \ # printf pattern to input frames
-f 1 -l 20 \                         # first and last frames
-sigma 20 \                          # noise std. dev.
-px1 4 -pt1 4 \                      # spatio-temporal patch size, step 1
-px2 4 -pt2 4 \                      # spatio-temporal patch size, step 2
-np1 1 -np2 1 \                      # number of similar patches, steps 1 & 2
-b1 3 \                              # noise correction factor beta, step 1
-deno output/path/final_%03d.png \   # [optional] where to store final estimate
-bsic output/path/basic_%03d.png \   # [optional] where to store basic estimate
-diff output/path/diffe_%03d.png     # [optional] where to store error sequence 

To denoise images

bin/vnlbayes \
-i puto/el/que/lee.png \ # printf pattern to input image
-sigma 20 \              # noise std. dev.
-px1 8 \                 # patch size, step 1
-px2 8 \                 # patch size, step 2
-np1 1 -np2 1 \          # number of similar patches, steps 1 & 2
-b1 3 \                  # noise correction factor beta, step 1
-deno output_final.png \ # [optional] where to store final estimate
-bsic output_basic.png \ # [optional] where to store basic estimate
-diff output_diffe.png   # [optional] where to store error image 

NOTE: The noise correction factor multiplies the noise std. dev. used in 
the MAP computation. This can be used to control the denoising strength.
For example, when using hard-thresholding the threshold is set to 
(b1*b1*sigma*sigma). 

Type bin/vnlbayes --help for more usage instructions. 

4/
A number of variants and combinations can be controlled by
defining/undefining a set of compilation flags (sorry) available in 

             src/NlBayes/VideoNLBayes.cpp

Please contact the authors for more information about them.

===============
ABOUT THIS FILE
===============

Copyright 2015 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

