/*
 * Copyright (c) 2016, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef LIBDCT_H_INCLUDED
#define LIBDCT_H_INCLUDED

#include <vector>
#include <fftw3.h>

struct DCTThreadsHandler
{
	// workspaces for DCT transforms in each thread
	fftwf_plan plan_forward [100];
	fftwf_plan plan_backward[100];
	float *dataspace[100];
	float *datafreq [100];

	// DCT bases for matrix product implementation 
	std::vector<float> basis_t;
	std::vector<float> basis_y;
	std::vector<float> basis_x;

	std::vector<float> dataspace_v[100];
	std::vector<float> datafreq_v [100];

	// size of the signals
	int width;
	int height;
	int frames;
	int nsignals;
	int nthreads;

	// empty constructors/destructors because construction
	// desctruction has to be handled manually
	DCTThreadsHandler(void) : nsignals(0) {};
	~DCTThreadsHandler(void) {};

	void init(int w, int h, int f, int n, int nthreads = 1);
	void destroy(void);

	// apply
	void forward(std::vector<float>& p);
	void inverse(std::vector<float>& p);

	private:
		enum MethodType {FFTW, MATPROD};
		MethodType method;
};

#endif //LIBDCT_H_INCLUDED
