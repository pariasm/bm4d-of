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

#include "dct_handler.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "matrix_funs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void DCTThreadsHandler::init(
	int w
,	int h
,	int f
,	int n
,	int t
){
	width = w;
	height = h;
	frames = f;
	nsignals = n;
	nthreads = t;

#ifdef _OPENMP 
	nthreads = std::min(nthreads, omp_get_max_threads());
	if (nthreads > 100)
	{
		fprintf(stderr,"Error: DCTThreadsHandler is hard-coded for a maximum of 100 threads\n");
		exit(1);
	}
#else
	if (nthreads > 1)
	{
		fprintf(stderr,"Error: DCTThreadsHandler can't handle %d threads (no OpenMP)\n", nthreads);
		exit(1);
	}
#endif


   unsigned int N = width * height * frames * nsignals;

	// define method based on patch size
//	method = (width * height * frames < 32) ? MATPROD : FFTW;
	method = FFTW;
//	method = MATPROD; // FIXME: MATPROD IS NOT WORKING!

//	fprintf(stderr, "init DCT for %d thread - %d x %d x %d ~ %d\n",nthreads, width, height, frames, nsignals);

	switch (method)
	{
		case FFTW:

//			fprintf(stderr, "init fftw for large patches\n");

			for (int i = 0; i < nthreads; i++)
			{
				dataspace[i] = (float*)fftwf_malloc(sizeof(float) * N);
				datafreq [i] = (float*)fftwf_malloc(sizeof(float) * N);
		
				int sz[] = {frames, height, width};
		
				fftwf_r2r_kind dct[] = {FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10};
				plan_forward[i] = fftwf_plan_many_r2r(3, sz, nsignals,
						dataspace[i], NULL, 1, width * height * frames, 
						datafreq [i], NULL, 1, width * height * frames,
//						dataspace[i], NULL, nsignals, 1, 
//						datafreq [i], NULL, nsignals, 1,
						dct, FFTW_ESTIMATE);
		
				fftwf_r2r_kind idct[] = {FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01};
				plan_backward[i] = fftwf_plan_many_r2r(3, sz, nsignals,
						datafreq [i], NULL, 1, width * height * frames,
						dataspace[i], NULL, 1, width * height * frames,
//						datafreq [i], NULL, nsignals, 1,
//						dataspace[i], NULL, nsignals, 1,
						idct, FFTW_ESTIMATE);
			}
			break;
			
		case MATPROD:

//			fprintf(stderr, "init matprod dct for small patches\n");

			for (int i = 0; i < nthreads; i++)
			{
				dataspace_v[i].resize(N);
				datafreq_v [i].resize(N);
			}

			// 1D DCT basis for signals of length width
			basis_x.resize(width  * width );
			for (unsigned kx = 0, i = 0; kx < width; kx++)
			for (unsigned nx = 0       ; nx < width; nx++, i++)
				basis_x[i] = cos( (0.5 + (float)nx)*  (float)kx * M_PI / (float)width )
				           * sqrt( (kx == 0 ? 1. : 2.) / (float)width ); 

			// 1D DCT basis for signals of length height
			basis_y.resize(height * height);
			for (unsigned ky = 0, i = 0; ky < height; ky++)
			for (unsigned ny = 0       ; ny < height; ny++, i++)
				basis_y[i] = cos( (0.5 + (float)ny)*  (float)ky * M_PI / (float)height )
				           * sqrt( (ky == 0 ? 1. : 2.) / (float)height ); 

			// 1D DCT basis for signals of length frames
			basis_t.resize(frames * frames);
			for (unsigned kt = 0, i = 0; kt < frames; kt++)
			for (unsigned nt = 0       ; nt < frames; nt++, i++)
//				basis_t[i] = nt == kt ? 1 : 0; 
				basis_t[i] = cos( (0.5 + (float)nt)*  (float)kt * M_PI / (float)frames )
				           * sqrt( (kt == 0 ? 1. : 2.) / (float)frames ); 

			break;
	}

}

void DCTThreadsHandler::destroy(void)
{
	if (nsignals)
	{
		for (int i = 0; i < nthreads; i++)
		{
			fftwf_free(dataspace[i]);
			fftwf_free(datafreq [i]);
			fftwf_destroy_plan(plan_forward [i]);
			fftwf_destroy_plan(plan_backward[i]);
		}
	}
}

void DCTThreadsHandler::forward(
	std::vector<float>& patch
){
   int tid = 0;
#ifdef _OPENMP 
   tid = omp_get_thread_num();
#endif

	int N = width * height * frames * nsignals;

	if (N == 0)
		fprintf(stderr, "Attempting to use a uninitialized DCTThreadsHandler.\n");

	switch (method)
	{
		/*
		case MATPROD: // compute dct via separable matrix products
		{
			// we assume a row-major memory layout for the patch
			const int whf = width * height * frames;
			const int wh  = width * height         ;
			const int wf  = width *          frames;
			
			// transform in x
			productMatrix(datafreq_v[tid], basis_x, patch,
							  width, nsignals*height*frames, width,
							  true, false);

			if (height == 1) // copy to patch and quit (assumes frames == 1)
			{
				patch = datafreq_v[tid];
				break;
			}

			// transform in y
			std::vector<float>::iterator p = datafreq_v[tid].begin();
			for (int n = 0; n < nsignals; ++n)
			{
				float *d = dataspace_v[tid].data() + n*whf;
				for (int t = 0; t < frames; ++t)
				for (int y = 0; y < height; ++y)
				for (int x = 0; x < width ; ++x)
					d[t*wh + x*height + y] = *p++;
			}

			productMatrix(datafreq_v[tid], basis_y, dataspace_v[tid],
							  height, nsignals*width*frames, height,
							  true, false);
			
			if (frames == 1) // copy to patch and quit
			{
				float *d = datafreq_v[tid].data();
				for (int n = 0; n < nsignals; ++n)
				{
					float *p = patch.data() + n*whf;
					for (int t = 0; t < frames; ++t)
					for (int x = 0; x < width ; ++x)
					for (int y = 0; y < height; ++y)
						p[t*wh + y*width + x] = *d++;
				}

				break;
			}

			// transform in t
			p = datafreq_v[tid].begin();
			for (int n = 0; n < nsignals; ++n)
			{
				float *d = dataspace_v[tid].data() + n*whf;
				for (int t = 0; t < frames; ++t)
				for (int x = 0; x < width ; ++x)
				for (int y = 0; y < height; ++y)
					d[y*wf + x*frames + t] = *p++;
			}

			productMatrix(datafreq_v[tid], basis_t, dataspace_v[tid],
							  frames, nsignals*width*height, frames,
							  true, false);
		
			// copy to patch and quit
			float *d = datafreq_v[tid].data();
			for (int n = 0; n < nsignals; ++n)
			{
				float *p = patch.data() + n*whf;
				for (int y = 0; y < height; ++y)
				for (int x = 0; x < width ; ++x)
				for (int t = 0; t < frames; ++t)
					p[t*wh + y*width + x] = *d++;
			}
		
			break;
		}
		*/

		case MATPROD: // compute dct via separable matrix products
		{
			// we assume a row-major memory layout for the patch
			const int whf = width * height * frames;
			const int wh  = width * height         ;
			const int wf  = width *          frames;
			
			// transform in x
			productMatrix(datafreq_v[tid], basis_x, patch,
							  width, nsignals*height*frames, width,
							  true, false);

			if (height == 1) // copy to patch and quit (assumes frames == 1)
			{
				patch = datafreq_v[tid];
				break;
			}

			// transform in y
			const int nf = frames*nsignals;
			for (int i = 0; i < nf; ++i)
			for (int yo = 0; yo < height; ++yo)
			{
				float *pstart = patch.data() + i*wh + yo*width;
				float *pstop  = pstart + width;
				for (int yb = 0; yb < height; ++yb)
				{
					float *d = datafreq_v[tid].data() + i*wh + yb*width;
					const float b = basis_y[yo*height + yb];
					for (float *p = pstart; p < pstop ; ++p, ++d)
						*p += *d * b; 
				} 
			}

//			const int nf = frames*nsignals;
//			for (int i = 0; i < nf; ++i)
//			for (int yo = 0; yo < height; ++yo)
//			{
//				float *pstart = patch.data() + i*wh + yo*width;
//				float *pstop  = pstart + width;
//
//				float *bstart = basis_y.data() + yo*height;
//				float *bstop  = bstart + height;
//
//				float *di = datafreq_v[tid].data() + i*wh;
//
//				for (const float *b = 0; b < bstop; ++b, di += width)
//				{
//					float *d = di;
//					for (float *p = pstart; p < pstop ; ++p, ++d)
//						*p += *d * *b; 
//				} 
//			}

			if (frames == 1) break;

			// transform in t
			std::vector<float>::iterator p = datafreq_v[tid].begin();
			for (int n = 0; n < nsignals; ++n)
			{
				float *d = dataspace_v[tid].data() + n*whf;
				for (int t = 0; t < frames; ++t)
				for (int x = 0; x < width ; ++x)
				for (int y = 0; y < height; ++y)
					d[y*wf + x*frames + t] = *p++;
			}

			productMatrix(datafreq_v[tid], basis_t, dataspace_v[tid],
							  frames, nsignals*width*height, frames,
							  true, false);
		
			// copy to patch and quit
			float *d = datafreq_v[tid].data();
			for (int n = 0; n < nsignals; ++n)
			{
				float *p = patch.data() + n*whf;
				for (int y = 0; y < height; ++y)
				for (int x = 0; x < width ; ++x)
				for (int t = 0; t < frames; ++t)
					p[t*wh + y*width + x] = *d++;
			}
		
			break;
		}

		case FFTW: // compute dct via fftw
		{
			// copy and compute unnormalized dct
			for (int i = 0; i < N; i++) dataspace[tid][i] = patch[i];

			fftwf_execute(plan_forward[tid]);

			// copy and orthonormalize
			float norm   = 1.0/sqrt(8.0*float(width*height*frames));
			float isqrt2 = 1.f/sqrt(2.0);

			for (int i = 0; i < N; i++) patch[i] = datafreq[tid][i] * norm;

			int whf = width * height * frames;
			int wh  = width * height;
			int w   = width;
			for (int n = 0; n < nsignals; n++)
			{
				for (int t = 0; t < frames; t++)
				for (int y = 0; y < height; y++)
					patch[n*whf + t*wh + y*w] *= isqrt2;

				for (int t = 0; t < frames; t++)
				for (int x = 0; x < width ; x++)
					patch[n*whf + t*wh + x] *= isqrt2;

				for (int y = 0; y < height; y++)
				for (int x = 0; x < width ; x++)
					patch[n*whf + y*w + x] *= isqrt2;
			}

			break;
		}
	}
}

void DCTThreadsHandler::inverse(
	std::vector<float>& patch
){
   int tid = 0;
#ifdef _OPENMP 
   tid = omp_get_thread_num();
#endif

	int N = width * height * frames * nsignals;

	if (N == 0)
		fprintf(stderr, "Attempting to use a uninitialized DCTThreadsHandler.\n");

	switch (method)
	{
		case MATPROD: // compute dct via separable matrix products
		{
			// we assume a row-major memory layout for the patch
			const int whf = width * height * frames;
			const int wh  = width * height         ;
			const int wf  = width *          frames;
			
			// transform in x
			productMatrix(datafreq_v[tid], basis_x, patch,
							  width, nsignals*height*frames, width,
							  false, false);

			if (height == 1) // copy to patch and quit (assumes frames == 1)
			{
				patch = datafreq_v[tid];
				break;
			}

			// transform in y
			std::vector<float>::iterator p = datafreq_v[tid].begin();
			for (int n = 0; n < nsignals; ++n)
			{
				float *d = dataspace_v[tid].data() + n*whf;
				for (int t = 0; t < frames; ++t)
				for (int y = 0; y < height; ++y)
				for (int x = 0; x < width ; ++x)
					d[t*wh + x*height + y] = *p++;
			}

			productMatrix(datafreq_v[tid], basis_y, dataspace_v[tid],
							  height, nsignals*width*frames, height,
							  false, false);
			
			if (frames == 1) // copy to patch and quit
			{
				float *d = datafreq_v[tid].data();
				for (int n = 0; n < nsignals; ++n)
				{
					float *p = patch.data() + n*whf;
					for (int t = 0; t < frames; ++t)
					for (int x = 0; x < width ; ++x)
					for (int y = 0; y < height; ++y)
						p[t*wh + y*width + x] = *d++;
				}

				break;
			}

			// transform in t
			p = datafreq_v[tid].begin();
			for (int n = 0; n < nsignals; ++n)
			{
				float *d = dataspace_v[tid].data() + n*whf;
				for (int t = 0; t < frames; ++t)
				for (int x = 0; x < width ; ++x)
				for (int y = 0; y < height; ++y)
					d[y*wf + x*frames + t] = *p++;
			}

			productMatrix(datafreq_v[tid], basis_t, dataspace_v[tid],
							  frames, nsignals*width*height, frames,
							  false, false);
		
			// copy to patch and quit
			float *d = datafreq_v[tid].data();
			for (int n = 0; n < nsignals; ++n)
			{
				float *p = patch.data() + n*whf;
				for (int y = 0; y < height; ++y)
				for (int x = 0; x < width ; ++x)
				for (int t = 0; t < frames; ++t)
					p[t*wh + y*width + x] = *d++;
			}
		
			break;
		}

		case FFTW: // compute dct via fftw
		{
			// normalize
			float norm  = 1.0/sqrt(8.0*float(width*height*frames));
			float sqrt2 = sqrt(2.0);

			int whf = width * height * frames;
			int wh  = width * height;
			int w   = width;
			for (int n = 0; n < nsignals; n++)
			{
				for (int t = 0; t < frames; t++)
				for (int y = 0; y < height; y++)
					patch[n*whf + t*wh + y*w] *= sqrt2;

				for (int t = 0; t < frames; t++)
				for (int x = 0; x < width ; x++)
					patch[n*whf + t*wh + x] *= sqrt2;

				for (int y = 0; y < height; y++)
				for (int x = 0; x < width ; x++)
					patch[n*whf + y*w + x] *= sqrt2;
			}

			for (int i = 0; i < N; i++) datafreq[tid][i] = patch[i] * norm;

			fftwf_execute(plan_backward[tid]);

			for (int i=0; i < N; i++) patch[i] = dataspace[tid][i];
		}
	}
}


