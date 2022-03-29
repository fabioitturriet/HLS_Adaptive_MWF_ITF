#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <ccomplex>
#include <stdio.h>


#include "fft_top.h"

#define NUM_CHAN 6

void mainFilter(float gama, float beta, int & ct_frame, std::complex<float> wm, float weighting[256], int estimation_mode, float frame_in_tmp_in[256][NUM_CHAN], int frame_vad[256], float & conty, float & contv, std::complex<float> (&rv)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&ry)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&rx)[NUM_CHAN][NUM_CHAN][256], float qL[NUM_CHAN], float qR[NUM_CHAN], std::complex<float> (&frame_out_tmp_in_aux)[256][2], std::complex<float> (&w)[2*NUM_CHAN][256]);
void mainFilterForResults(int ct_frame, std::complex<float> wm, float weighting[256], float frame_in_tmp_in[256][NUM_CHAN], float qL[NUM_CHAN], float qR[NUM_CHAN], std::complex<float> (&frame_out_frq_out_puro)[256][2], std::complex<float> (&frame_out_frq_out_filtrado)[256][2], std::complex<float> w[12][256]);

void mwf_itf(float gama, float beta, std::complex<float> initial_coeffs[12][256], std::complex<float> rv[NUM_CHAN][NUM_CHAN][256], std::complex<float> ry[NUM_CHAN][NUM_CHAN][256], std::complex<float> rx[NUM_CHAN][NUM_CHAN][256], int m, float qL[NUM_CHAN], float qR[NUM_CHAN]);
void corr_matrix_estimation(int estimation_mode, int length_fft, int channel, std::complex<float> frame_y[256][NUM_CHAN], std::complex<float> (&rv)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&ry)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&rx)[NUM_CHAN][NUM_CHAN][256], int vad, float & conty, float & contv, float smooth_coef);
