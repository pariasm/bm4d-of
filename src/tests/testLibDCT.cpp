/*
 * Copyright (c) 2015, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cassert>

#include <string>
#include <sstream>

#include "Utilities.h"
#include "LibDCT.h"
#include "LibVideoT.hpp"
#include "LibMatrix.h"


/**
 * @file   testLibVideoT.cpp
 * @brief  Executable file to test libVideoT
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

int main(int argc, char **argv)
{
	//! Check if there is the right call for the algorithm
	if (argc < 4)
	{
		fprintf(stdout, "Usage: %s direction in-path out-path"
		                "first-frame last-frame\n", argv[0]);
		return EXIT_FAILURE;
	}

	//! Get command line inputs
	const int direction = atoi(argv[1]);
	const char*  input_path =  argv[2] ;
	const char* output_path =  argv[3] ;
	const int first_frame = (argc > 4) ? atoi(argv[4]) : 0;
	const int  last_frame = (argc > 5) ? atoi(argv[5]) : first_frame;
	const int frame_step = 1;

	if ((direction != 1) && (direction != -1))
	{
		fprintf(stdout, "DCT's direction must be 1 or -1\n");
		return EXIT_FAILURE;
	}

	//! Load video
	Video<float> v(input_path, first_frame, last_frame, frame_step);

	printMatrix(v.data,v.sz.height, v.sz.width, "input.asc");

	//! DCT
	DCTThreadsHandler dct;

	dct.init(v.sz.width, v.sz.height, v.sz.frames, v.sz.channels);
	if (direction == 1) dct.forward(v.data);
	else                dct.inverse(v.data);
	dct.destroy();

	printMatrix(v.data,v.sz.height, v.sz.width, "output.asc");

	//! Save video
	v.saveVideo(output_path, first_frame, frame_step);

	return EXIT_SUCCESS;
}
