/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <complex>
#include <gnuradio/gr_complex.h>
#include <fftw3.h>
#include "framing.h"
#include <stdlib.h>

// OFDM configurations
#define NUM_SUBCARRIERS     64
#define USE_FFTW            true

#define LFSR_LENGTH   12
#define LFSR_GEN_POLY 010123

int main(int argc, char **argv)
{
  unsigned int M = NUM_SUBCARRIERS;
  msequence ms = msequence_create(LFSR_LENGTH,
				  LFSR_GEN_POLY,
				  1);
  gr_complex * x;
  gr_complex * X;
  gr_complex * Y;
  unsigned char * c;
  float g;
  unsigned int M_S1 = 0;
  fftwf_plan fft, ifft;
  x = (gr_complex *)fftwf_malloc(sizeof(gr_complex)*M);
  X = (gr_complex *)fftwf_malloc(sizeof(gr_complex)*M);
  Y = (gr_complex *)fftwf_malloc(sizeof(gr_complex)*M);
  c = (unsigned char *)malloc(sizeof(unsigned char)*M);
  #if USE_FFTW
    ifft = fftwf_plan_dft_1d(M,
  	  		     reinterpret_cast<fftwf_complex *>(X),
			     reinterpret_cast<fftwf_complex *>(x),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fft  = fftwf_plan_dft_1d(M,
			     reinterpret_cast<fftwf_complex *>(x),
			     reinterpret_cast<fftwf_complex *>(Y),
                             FFTW_FORWARD,
                             FFTW_ESTIMATE);
  #endif

  ofdmframe_init_default_sctype(c, M);
  printf("******* SUBCARRIER ALLOCATION ********\n[");
  for(unsigned int i = 0; i < M; i++) {
    switch(c[i]) {
      case OFDMFRAME_SCTYPE_NULL:     printf(".");    break;
      case OFDMFRAME_SCTYPE_PILOT:    printf("|");    break;
      case OFDMFRAME_SCTYPE_DATA:     printf("+");    break;
    }
  }
  printf("]\n******* SUBCARRIER ALLOCATION ********\n");
  ofdmframe_print_sctype(c, M);
  for(unsigned int i = 0; i < M; i++) {
    if(c[i] == OFDMFRAME_SCTYPE_NULL)
      X[i] = 0.0f;
    else {
      X[i] = ((msequence_generate_symbol(ms, 3) & 0x01) ? 1.0f : -1.0f);
      M_S1++;
    }
  }
  g = 1.0f / M;

  #if USE_FFTW
    fftwf_execute(ifft);
    for(unsigned int i = 0; i < M; i++)
      x[i] *= g;
    fftwf_execute(fft);
  #else
    fft_run(M, X, x, LIQUID_FFT_BACKWARD, 0);
    for(unsigned int i = 0; i < M; i++)
      x[i] *= g;
    fft_run(M, x, Y, LIQUID_FFT_FORWARD, 0);
  #endif

  printf("      X\t\t\t      Y\n");
  for(unsigned int i = 0; i < M; i++) {
    printf("%3.3f\t%3.3fi\t\t%3.3f\t%3.3f\n",
	   std::real(X[i]),
	   std::imag(X[i]),
	   std::real(Y[i]),
	   std::imag(Y[i]));
  }

  #if USE_FFTW
    fftwf_destroy_plan(fft);
    fftwf_destroy_plan(ifft);
  #endif
  free(c);
  fftwf_free(x);
  fftwf_free(X);
  fftwf_free(Y);
  return 0;
}
