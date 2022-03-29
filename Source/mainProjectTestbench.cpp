/*******************************************************************************
Vendor: Xilinx 
Associated Filename: fft_tb.cpp
Purpose: Xilinx FFT IP-XACT IP in Vivado HLS
Revision History: September 26, 2013 - initial release
                                                
*******************************************************************************
#-  (c) Copyright 2011-2018 Xilinx, Inc. All rights reserved.
#-
#-  This file contains confidential and proprietary information
#-  of Xilinx, Inc. and is protected under U.S. and
#-  international copyright and other intellectual property
#-  laws.
#-
#-  DISCLAIMER
#-  This disclaimer is not a license and does not grant any
#-  rights to the materials distributed herewith. Except as
#-  otherwise provided in a valid license issued to you by
#-  Xilinx, and to the maximum extent permitted by applicable
#-  law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND
#-  WITH ALL FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES
#-  AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, INCLUDING
#-  BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-
#-  INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE; and
#-  (2) Xilinx shall not be liable (whether in contract or tort,
#-  including negligence, or under any other theory of
#-  liability) for any loss or damage of any kind or nature
#-  related to, arising under or in connection with these
#-  materials, including for any direct, or any indirect,
#-  special, incidental, or consequential loss or damage
#-  (including loss of data, profits, goodwill, or any type of
#-  loss or damage suffered as a result of any action brought
#-  by a third party) even if such damage or loss was
#-  reasonably foreseeable or Xilinx had been advised of the
#-  possibility of the same.
#-
#-  CRITICAL APPLICATIONS
#-  Xilinx products are not designed or intended to be fail-
#-  safe, or for use in any application requiring fail-safe
#-  performance, such as life-support or safety devices or
#-  systems, Class III medical devices, nuclear facilities,
#-  applications related to the deployment of airbags, or any
#-  other applications that could lead to death, personal
#-  injury, or severe property or environmental damage
#-  (individually and collectively, "Critical
#-  Applications"). Customer assumes the sole risk and
#-  liability of any use of Xilinx products in Critical
#-  Applications, subject only to applicable laws and
#-  regulations governing limitations on product liability.
#-
#-  THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS
#-  PART OF THIS FILE AT ALL TIMES. 
#- ************************************************************************


This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/

# define M_PI           3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
//#include "fft_top.h"
#include "mainProject.hpp"
#include <stdio.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <fstream>
#include <string>
#include <sstream>


using namespace std;

#define NUM_CHAN 6    //Total Number of Microphones
#define BUF_SIZE 64

float sciToDub(const string& str) {

   stringstream ss(str);
   float d = 0;
   ss >> d;

   if (ss.fail()) {
      string s = "Unable to format ";
      s += str;
      s += " as a number!";
      throw (s);
   }

   return (d);
}

bool compare_float(float x, float y, float epsilon = 0.001f){
   if(fabs(x - y) < epsilon)
      return true; //they are same
      return false; //they are not same
}

int main()
{

	std::cout.precision(20);
	int tamSign = 188879; // tamanho sinal

	std::string line, token;
	std::ifstream infilePeso("peso.txt");


	ifstream* inputAudio = new ifstream;
	ifstream* inputAudioVAD = new ifstream;

	int countline = 0;
//	int NUM_CHAN = 6;

	countline = 0;

	float weighting[256];
	while (std::getline(infilePeso, line)){
		std::istringstream iss(line);
		if(countline>255)
			exit(0);

		std::getline(iss, token, '\n');
		float numberDub = sciToDub(token);
		weighting[countline] = numberDub;
		countline++;
	}
	float gama[2] = {0.0001, 0.0006};
	float beta = 0.015;
	std::complex<float> wm(0.999698818696204,0.024541228522912);

	for (int numTeste = 0; numTeste < 1; numTeste++){
		for (int countGama = 0; countGama < 1; countGama++){

			if(numTeste == 0)
			{
				inputAudio->open("TestBenchF+60.txt");
				inputAudioVAD->open("VADTestBenchF+60.txt");
			}
			else
			{
				inputAudio->open("TestBenchF-60.txt");
				inputAudioVAD->open("VADTestBenchF-60.txt");
			}
		// INICIALIZACAO

			std::complex<float> rv[6][6][256], ry[6][6][256], rx[6][6][256];
			int m=256;
			float qL[6], qR[6];
			for (int i=0; i<NUM_CHAN;i++)
			{
				qL[i] = 0;
				qR[i] = 0;
			}

			qL[0] = 1;
			qR[((NUM_CHAN/2)+1)-1] = 1;
			int length_fft = 256;
			int channel = 6;
			int sumVad = 0;
			std::complex<float> frame_y[256][6];
			int frameEndPriorEst = 200;
			int j = 0;

			for(int bin = 0; bin < NUM_CHAN; bin++)
			{
				for(int i = 0; i < NUM_CHAN; ++i)
				{
					for(int j = 0; j < 256; ++j)
					{
						std::complex<float> entry(0,0);
						rv[bin][i][j] = entry;
						rx[bin][i][j] = entry;
						ry[bin][i][j] = entry;
					}
				}
			}
			// FIM INICIALIZACAO


			// IMPORTANTES
			float input[6];
			float frame_in_tmp_in[256][6];
			//

			std::complex<float> resultFFT[256];

			float realPart=0.0, imagPart=0.0;

			// FIM MEU CODIGO

			int ct_frame = 0;
			int frame_vad[256];
			//frame_vad[0] = 1;
			float conty = 0;
			float contv = 0;
			int estimation_mode = 0;
			std::complex<float> frame_out_tmp_in[256][2];
			std::complex<float> frame_out_tmp_in_aux[256][2];
			std::complex<float> w[6*2][256];

			for (int i = 0; i < 256; i++)
			{
				 w[0][i] = 1;
				 w[9][i] = 1; // era 9 para 6 mics
			}

			ofstream audioSaida;
			int countVAD = 0;
			string saidaNome = "audioSaida_Beta" + std::to_string(beta) + "_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			audioSaida.open (saidaNome);

			ofstream wSaida;
			string wSaidaNome = "saidaW_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			wSaida.open (wSaidaNome);

			for (int sampleCount = 0; sampleCount < tamSign; sampleCount++)
			{
				for (int k =0; k < 255; k++)
					for (int g =0; g < NUM_CHAN; g++)
						frame_in_tmp_in[k][g] = frame_in_tmp_in[k+1][g];

				std::getline(*inputAudio, line);
				std::istringstream iss(line);

				for (int i = 0; i < NUM_CHAN; i++) {
					std::getline(iss, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}


				for(int i = 255; i > 0; i--) {
					frame_vad[i] = frame_vad[i-1];
				}

				std::getline(*inputAudioVAD, line);
				std::istringstream iss2(line);
				std::getline(iss2, token, '\n');
				int numberDubInt = (int)sciToDub(token);
				frame_vad[0] = numberDubInt;


				for (int g = 0; g < NUM_CHAN; g++)
					frame_in_tmp_in[255][g] = input[g];

				for (int k =0; k < 255; k++)
					for (int g =0; g < 2; g++)
						frame_out_tmp_in[k][g] = frame_out_tmp_in[k+1][g];

				for (int g =0; g < 2; g++)
					frame_out_tmp_in[255][g] = 0;

				if (((sampleCount+1) % 64) == 0)
				{
					mainFilter(gama[countGama], beta, ct_frame, wm, weighting, estimation_mode, frame_in_tmp_in, frame_vad, conty, contv, rv, ry, rx, qL, qR, frame_out_tmp_in_aux, w);

					for (int i = 0; i < 256; i++)
						for (int k = 0; k < 2; k++)
						{
							frame_out_tmp_in[i][k] =frame_out_tmp_in[i][k] + frame_out_tmp_in_aux[i][k];
						}

					int x = 0;
					for (int outCount = 0; outCount < 64; outCount++)
					{
						for (int j = 0; j < 2; j++)
						{
							audioSaida << std::fixed << std::setprecision(16) << frame_out_tmp_in[outCount][j].real();
							if (j==0)
								audioSaida << "\t";
						}
						audioSaida << "\n";
						x++;
					}
					for (int wCount = 0; wCount < 256; wCount++)
					{
						for (int wCountMic = 0; wCountMic < 12; wCountMic++)
						{
							wSaida <<  std::fixed << std::setprecision(16) << w[wCountMic][wCount].real() << " " << w[wCountMic][wCount].imag() << "\t";
						}
						wSaida << "\n";
					}
				}
			}
			wSaida.close();
		    audioSaida.close();
			inputAudio->close();
			inputAudioVAD->close();
		}
	}
	/*
    int error = 0;
    int errorGeral = 0;
    int error4Casa = 0;
    int error3Casa = 0;
    int error2Casa = 0;
    int errorAcima = 0;

    for (int numTeste = 0; numTeste < 2; numTeste++){
		for (int countGama = 0; countGama < 1; countGama++){
			std::ifstream inW("saidaW_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt");
			std::ifstream inFala("fala.txt");
			std::ifstream inRuido;

			if (numTeste == 0)
				inRuido.open("ruidoF+60.txt");
			else
				inRuido.open("ruidoF-60.txt");



			ofstream snrPuro;
			string snrPuroNome = "snrPuro_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			snrPuro.open (snrPuroNome);
			ofstream snrFiltrado;
			string snrFiltradoNome = "snrFiltrado_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			snrFiltrado.open (snrFiltradoNome);
			ofstream deltaSnr;
			string deltaSnrNome = "deltaSnr_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			deltaSnr.open (deltaSnrNome);
			ofstream ILD_lf;
			string ILD_lfNome = "ild_lf_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			ILD_lf.open (ILD_lfNome);
			ofstream ild;
			string ildNome = "ild_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			ild.open (ildNome);
			ofstream deltaILD;
			string deltaILDNome = "delta_ILD_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			deltaILD.open (deltaILDNome);
			ofstream deltaITD;
			string deltaITDNome = "delta_ITD_Gama" + std::to_string(gama[countGama]) + "_NumTeste" + std::to_string(numTeste) + ".txt";
			deltaITD.open (deltaITDNome);


			// IMPORTANTES
			float input[6];
			float frame_in_tmp_inFala[256][6];
			float frame_in_tmp_inRuido[256][6];
			for (int i=0; i<256;i++)
				for (int j=0; j<6;j++)
				{
					frame_in_tmp_inFala[i][j] = 0;
					frame_in_tmp_inRuido[i][j] = 0;
				}

			int ct_frame = 0;
			float qL[6], qR[6];
			for (int i=0; i<5;i++)
			{
				qL[i] = 0;
				qR[i] = 0;
			}

			qL[0] = 1;
			qR[((6/2)+1)-1] = 1;
			std::complex<float> frame_out_frq_out_puroFala[256][2];
			std::complex<float> frame_out_frq_out_filtradoFala[256][2];
			std::complex<float> frame_out_frq_out_puroRuido[256][2];
			std::complex<float> frame_out_frq_out_filtradoRuido[256][2];
			std::complex<float> wSaida[12][256];

			float snrPuroL = 0;
			float snrPuroR = 0;
			float snrFiltradoL = 0;
			float snrFiltradoR = 0;

			float aux_sp = 0;
			float aux_no = 0;

			float vec_in = 0;
			float vec_out = 0;

			float ILD_in_LF = 0;
			float ILD_out_LF = 0;

			float ild_sp = 0;
			float ild_no = 0;

			float Delta_ild_sp = 0;
			float Delta_ild_no = 0;

			float itd_sp[256];
			float itd_no[256];
			std::complex<float> itd_sp_o[256];
			std::complex<float> itd_no_o[256];
			std::complex<float> itd_sp_i[256];
			std::complex<float> itd_no_i[256];
			float ild_sp_vec[256];
			float ild_no_vec[256];

			float temp_itd_sp[128];
			float temp_itd_no[128];
			float Delta_itd_sp = 0;
			float Delta_itd_no = 0;

	//		minfrq = floor(0*M/FAM+1);
	//		aux_sp = (pow_frq(:,3).*pow_frq(:,5))./(pow_frq(:,1).*pow_frq(:,7));
	//		aux_no = (pow_frq(:,4).*pow_frq(:,6))./(pow_frq(:,2).*pow_frq(:,8));

	//		vec_in=pow_frq(:,2)./pow_frq(:,6);
	//		vec_out=pow_frq(:,4)./pow_frq(:,8);
	//
	//		ILD_in_LF =  10*log10(vec_in(1:maxfrq));
	//		ILD_out_LF = 10*log10(vec_out(1:maxfrq));
	//
	//
	//		ild_sp = abs(10*log10(aux_sp(minfrq:M/2)));
	//		ild_no = abs(10*log10(aux_no(minfrq:M/2)));
	//
	//		Delta_ild_sp = (1 / (M/2+1-minfrq) ) * sum(ild_sp);
	//		Delta_ild_no = (1 / (M/2+1-minfrq) ) * sum(ild_no);

			for (int sampleCount = 0; sampleCount < tamSign; sampleCount++)
			{
				for (int k =0; k < 255; k++)
					for (int g =0; g < 6; g++)
						frame_in_tmp_inFala[k][g] = frame_in_tmp_inFala[k+1][g];

				std::getline(inFala, line);
				std::istringstream iss(line);
				for (int i = 0; i < 6; i++) {
					std::getline(iss, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}


				for (int g = 0; g < 6; g++)
					frame_in_tmp_inFala[255][g] = input[g];

				for (int k =0; k < 255; k++)
					for (int g =0; g < 6; g++)
						frame_in_tmp_inRuido[k][g] = frame_in_tmp_inRuido[k+1][g];

				std::getline(inRuido, line);
				std::istringstream iss2(line);
				for (int i = 0; i < 6; i++) {
					std::getline(iss2, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}


				for (int g = 0; g < 6; g++)
					frame_in_tmp_inRuido[255][g] = input[g];

				if (((sampleCount+1) % 64) == 0)
				{
					ct_frame++;


					float realPart, imagPart;
					for (int j = 0; j < 256; j++){
						std::getline(inW, line);
						std::istringstream iss(line);
						for (int i = 0; i < 12; i++){
							std::getline(iss, token, ' ');
							realPart = sciToDub(token);
							std::getline(iss, token, '\t');
							imagPart = sciToDub(token);
							std::complex<float> entry(realPart, imagPart);
							wSaida[i][j] = entry;
						}
					}

					float somaFala = 0;
					float somaFalaRuido = 0;
					for (int k = 0; k < 256; k++)
					{
						for (int y = 0; y < 6; y++)
						{
							somaFala = somaFala + std::abs(frame_in_tmp_inFala[k][y]);
							somaFalaRuido = somaFalaRuido + std::abs(frame_in_tmp_inRuido[k][y]);
						}
					}

					if (somaFala > 0.0)
					{
						mainFilterForResults(ct_frame, wm, weighting, frame_in_tmp_inFala, qL, qR, frame_out_frq_out_puroFala, frame_out_frq_out_filtradoFala, wSaida);

						if (somaFalaRuido > 0.0)
							mainFilterForResults(ct_frame, wm, weighting, frame_in_tmp_inRuido, qL, qR, frame_out_frq_out_puroRuido, frame_out_frq_out_filtradoRuido, wSaida);

						float powerFreqPuroFala[2] ={0.0,0.0};
						float powerFreqFiltradoFala[2] = {0.0,0.0};
						float powerFreqPuroRuido[2] ={0.0,0.0};
						float powerFreqFiltradoRuido[2] = {0.0,0.0};

						for (int lado = 0; lado < 2; lado++)
						{
							float somaBinPuroFala = 0;
							float somaBinFiltradoFala = 0;
							float somaBinPuroRuido = 0;
							float somaBinFiltradoRuido = 0;

							for (int bin = 0; bin < 256; bin++)
							{
								somaBinPuroFala = somaBinPuroFala + pow(std::abs(frame_out_frq_out_puroFala[bin][lado]),2);
								somaBinFiltradoFala = somaBinFiltradoFala + pow(std::abs(frame_out_frq_out_filtradoFala[bin][lado]),2);
								somaBinPuroRuido = somaBinPuroRuido + pow(std::abs(frame_out_frq_out_puroRuido[bin][lado]),2);
								somaBinFiltradoRuido = somaBinFiltradoRuido + pow(std::abs(frame_out_frq_out_filtradoRuido[bin][lado]),2);


							}
							powerFreqPuroFala[lado] = somaBinPuroFala;
							powerFreqFiltradoFala[lado] = somaBinFiltradoFala;
							powerFreqPuroRuido[lado] = somaBinPuroRuido;
							powerFreqFiltradoRuido[lado] = somaBinFiltradoRuido;
						}


						for (int bin = 0; bin < 256; bin++)
						{
							itd_sp[bin] =  pow(std::arg(frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_filtradoFala[bin][1]) * frame_out_frq_out_puroFala[bin][1] * std::conj(frame_out_frq_out_puroFala[bin][0])),2);
							itd_no[bin] =  pow(std::arg(frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_filtradoRuido[bin][1]) * frame_out_frq_out_puroRuido[bin][1] * std::conj(frame_out_frq_out_puroRuido[bin][0])),2);
							itd_sp_o[bin] = frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_filtradoFala[bin][1]);
							itd_sp_i[bin] = frame_out_frq_out_puroFala[bin][0] * std::conj(frame_out_frq_out_puroFala[bin][1]);
							itd_no_i[bin] = frame_out_frq_out_puroRuido[bin][0] * std::conj(frame_out_frq_out_puroRuido[bin][1]);
							itd_no_o[bin] = frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_filtradoRuido[bin][1]);

							ild_sp_vec[bin] = pow( (log10((abs(frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_puroFala[bin][1])))  / (abs( frame_out_frq_out_filtradoFala[bin][1] * std::conj(frame_out_frq_out_puroFala[bin][0]))) ) ) , 2);
							ild_no_vec[bin] = pow( (log10((abs(frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_puroRuido[bin][1])))  / (abs( frame_out_frq_out_filtradoRuido[bin][1] * std::conj(frame_out_frq_out_puroRuido[bin][0]))) ) ) , 2);

	//						ild_sp(bin) = ild_sp(bin) + ( log10(   abs(    out_frq(bin,3) * conj(out_frq(bin,5))    ) / abs(   out_frq(bin,7) * conj(out_frq(bin,1))  )    ) )^2;

	//						itd_sp(bin) = itd_sp(bin) + ( phase(    out_frq(bin,3) * conj(out_frq(bin,7)) * out_frq(bin,5) * conj(out_frq(bin,1))        ) )^2;

						}

						float soma_itd_no = 0;
						float soma_itd_sp = 0;
						for (int bin = 0; bin < 128; bin++)
						{
							temp_itd_sp[bin] = abs(std::arg(itd_sp_o[bin]) - std::arg(itd_sp_i[bin]));
							temp_itd_no[bin] = abs(std::arg(itd_no_o[bin]) - std::arg(itd_no_i[bin]));

							if (temp_itd_no[bin] > M_PI)
							{
								temp_itd_no[bin] = abs((2*M_PI)-temp_itd_no[bin]);
							}
							if (temp_itd_sp[bin] > M_PI)
							{
								temp_itd_sp[bin] = abs((2*M_PI)-temp_itd_sp[bin]);
							}

							soma_itd_no = soma_itd_no + temp_itd_no[bin];
							soma_itd_sp = soma_itd_sp + temp_itd_sp[bin];
						}


						Delta_itd_sp= (soma_itd_sp/128)/M_PI;
						Delta_itd_no= (soma_itd_no/128)/M_PI;


						if (powerFreqPuroFala[0] > 0.00001)
						{
							snrPuroL = 10 * log10(powerFreqPuroFala[0]/powerFreqPuroRuido[0]);
							snrPuroR = 10 * log10(powerFreqPuroFala[1]/powerFreqPuroRuido[1]);
							snrFiltradoL = 10 * log10(powerFreqFiltradoFala[0]/powerFreqFiltradoRuido[0]);
							snrFiltradoR = 10 * log10(powerFreqFiltradoFala[1]/powerFreqFiltradoRuido[1]);

							aux_sp = (powerFreqFiltradoFala[0]*powerFreqPuroFala[1]) / (powerFreqPuroFala[0] * powerFreqFiltradoFala[1]);
							aux_no = (powerFreqFiltradoRuido[0]*powerFreqPuroRuido[1]) / (powerFreqPuroRuido[0] * powerFreqFiltradoRuido[1]);;

							vec_in = powerFreqPuroRuido[0] / powerFreqPuroRuido[1];
							vec_out = powerFreqFiltradoRuido[0] / powerFreqFiltradoRuido[1];

							ILD_in_LF = 10*log10(vec_in); // ILD de entrada, ruido esquerdo / ruido direito
							ILD_out_LF = 10*log10(vec_out); // ILD de saida, ruido filtrado esquerdo / ruido filtrado direito

							ild_sp = abs(10*log10(aux_sp)); // ILD fala, ?
							ild_no = abs(10*log10(aux_no)); // ILD ruido, ?

							Delta_ild_sp = 0.33333333 * ild_sp; // Delta ILD fala, ?
							Delta_ild_no = 0.33333333 * ild_no; // Delta ILD ruido, ?
						}
						else
						{
							snrPuroL = -10;
							snrPuroR = -10;
							snrFiltradoL = -10;
							snrFiltradoR = -10;
						}

						//std::cout << "PURO= SNR L: " << snrPuroL << " SNR R: " << snrPuroR <<  "FILTRADO= SNR L: " << snrFiltradoL << " SNR R: " << snrFiltradoR << "\n";
					}
					else
					{
						snrPuroL = -10;
						snrPuroR = -10;
						snrFiltradoL = -10;
						snrFiltradoR = -10;

						//std::cout << "PURO= SNR L: " << snrPuroL << " SNR R: " << snrPuroR <<  "FILTRADO= SNR L: " << snrFiltradoL << " SNR R: " << snrFiltradoR << "\n";
					}

					snrPuro << std::fixed << std::setprecision(16) << snrPuroL << "\t" << snrPuroR << "\n";
					snrFiltrado << std::fixed << std::setprecision(16) << snrFiltradoL << "\t" << snrFiltradoR << "\n";
					deltaSnr << std::fixed << std::setprecision(16) << (snrFiltradoL-snrPuroL) << "\t" << (snrFiltradoR-snrPuroR) << "\n";
					ILD_lf << std::fixed << std::setprecision(16) << ILD_in_LF << "\t" << ILD_out_LF << "\n";
					ild << std::fixed << std::setprecision(16) << ild_sp << "\t" << ild_no << "\n";
					deltaILD << std::fixed << std::setprecision(16) << Delta_ild_sp << "\t" << Delta_ild_no << "\n";
					deltaITD << std::fixed << std::setprecision(16) << Delta_itd_sp << "\t" << Delta_itd_no << "\n";
				}
			}

			inW.close();
			inFala.close();
			inRuido.close();
			snrPuro.close();
			snrFiltrado.close();
			deltaSnr.close();
			ILD_lf.close();
			ild.close();
			deltaILD.close();
			deltaITD.close();
		}
    }
	 */

        return 0;
}

