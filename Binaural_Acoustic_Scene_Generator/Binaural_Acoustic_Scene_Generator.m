% -------------------------------------------------------------------------
% Binaural_Acoustic_Scene_Generator.m
%
% Generates binaural speech and voice activity information.
%
% [1] H. Kayser, S.D. Ewert, J. Anemüller, T. Rohdenburg, V. Hohmann, B. 
% Kollmeier (2009) Database of multichannel in-ear and behind-the-ear
% head-related and binaural room impulse responses. EURASIP Journal on
% Advances in Signal Processing, pages 1-10.
% http://medi.uni-oldenburg.de/hrir/      (v1.1)
%
% [2] IEEE Subcommittee (1969)IEEE Recommended practice for speech 
% quality measurements. IEEE Transactions on Audio and Electroacoustics
% AU-17(3):225-246. 
% 
% Inputs: Data from the acoustic scene
% Outputs: .wav and .txt files with microphones signals
% 
% Author    : Fábio Pires Itturriet & Diego Marques
% Version   : 1.0
% Created   : 27/09/2021 
% Last Modif:
% 
% -------------------------------------------------------------------------

% Cleaning previous workspace
clear all;
close all;
% Environment Definition ---------------------------------------------------
% Options = [Anechoic  Office_I  Office_II  Cafeteria  Courtyard]
% Befor select the enviroment, please read the documentation (/Doc)
ENVIRONMENT='Anechoic';       % Definition among five different scenarios

% parameters---------------------------------------------------------------
log   = 1;     % progress information: 1 - on, 2 - off
plots = 0;     % plots during the stages: 1 - on, 0 - off

% Position of the speech source -------------------------------------------
DTS =    300;   % speech distance  = [ 80 300 ]
ELS =      0;   % speech elevation = [ down -10 0 10 20 up ]
AZS =      0;   % speech azimuth   = { left -180 : 5 : 180 right }

% Position of the noise source---------------------------------------------
DTN =    300;    % noise distance   = [ 80 300 ]
ELN =      0;    % noise elevation  = [ -10 0 10 20 ]
AZN =    -60;    % noise azimuth    = [ -180 : 5 : 180 ]

% Speech files-------------------------------------------------------------
SpeechPFile_01='S_50_08_16kHz.wav';   % Male voice
SpeechPFile_02='ieee03f09_16kHz.wav'; % Female voice 

% Selected Speech File
SPEECH_FILE  = SpeechPFile_02;

% Noise files -------------------------------------------------------------
Noise_File_01= 'ICRA_No01_16kHz_12s.wav';
Noise_File_02= 'Cafeteria_Babble_16kHz_12s.wav';
Noise_File_03= 'Honda_Bike_Start_Up_2_16kHz_12s.wav';            
Noise_File_04= 'Street_Noise_Loizou_16kHz_12s.wav';
Noise_File_05= 'Cafeteria_Babble_Loizou_16kHz_12s.wav';             
Noise_File_06= 'FX_Bus_Engine_001_16kHz_12s.wav';
Noise_File_07= 'SpeechLike_A_16kHz_12s.wav';

% Selected Noise File
NOISE_FILE  = Noise_File_01;
          
% Signal-to-Noise Ratio (SNR) in dB ---------------------------------------
input_SNR=0;
SNS=0;  %

% Number of Microphones (L+R=Total)--------------------------------------------
M=2;         % Total of Binaural Hearing aids = [ 2(1-1) 4(2-2) 6(3-3) ]

% Normalization factor ----------------------------------------------------
FCT =   0.9;    % normalization factor (for wave recordings)

%**************************************************************************
%******************   Nothing to change from here   *********************** 
%**************************************************************************
% directories and file names  ---------------------------------------------
DIROUT    = './Results/';                      % output file directory
DIRSPEECH = './Data/Speech/';                  % input speech file directory
DIRHRIR   = './Data/HRIR_database_mat/';       % hrir file directory
DIRNOISE  = './Data/Noise/';                   % noise  file directory
FnamSpOut = strcat(DIROUT,'I1_SPEECH_Front_Mics','.wav'); 
FnamNoOut = strcat(DIROUT,'I2_NOISE_Front_Mics','.wav'); 
FnamInFront=strcat(DIROUT,'I3_INPUT_Front_Mics','.wav'); 

% speech  -----------------------------------------------------------------
[raw_speech,FAM] = audioread([DIRSPEECH,SPEECH_FILE]);

% noise  -----------------------------------------------------------------
[raw_noise,FAM_NO]  = audioread([DIRNOISE,NOISE_FILE]);

% Loading HRIRs [1] --------------------------------------------------------
old = cd;
cd(DIRHRIR)
hrir_speech = loadHRIR(ENVIRONMENT,DTS,ELS,AZS,'bte');  % odd channel is left side
hrir_noise  = loadHRIR(ENVIRONMENT,DTN,ELN,AZN,'bte');  % odd channel is left side
cd(old)
HRIRSP = resample(hrir_speech.data,FAM,48000);
HRIRNO = resample(hrir_noise.data,FAM,48000);
clear old hrir;

% Define the number of channel according to the number of Mics (M)
numchan = size(HRIRSP,2);

% length of convolved signals ---------------------------------------------
tamsign = length(raw_speech) + length(HRIRSP(:,1)) - 1;
tamsp_ori = tamsign + length(HRIRNO) - 1;

%**************************************************************************
% Fitting the selected Noise
%**************************************************************************
    
 tamno_raw = length(raw_noise);
 if ( (FAM_NO~=FAM) || (tamno_raw<tamsp_ori) )
      disp('ERROR - Different noise/speech sampling frequencies or lengths');
      pause;
 end       
 raw_noise = raw_noise(1:tamsp_ori); %Speech and noise with same size       
       
%**************************************************************************
% Convolution with hrir
%**************************************************************************

% convolved speech and noise ----------------------------------------------
speech_hri = zeros(tamsign,numchan);
noise_hri  = zeros(tamsign,numchan);

for cont = 1:numchan/2
    speech_hri( :, cont ) = conv( raw_speech, HRIRSP(:,2*(cont-1)+1), 'full' );
    speech_hri( :, cont+numchan/2 ) = conv( raw_speech, HRIRSP(:,2*cont), 'full' );
    noise_hri(:,cont) = conv(raw_noise,HRIRNO(:,2*(cont-1)+1),'valid');
    noise_hri(:,cont+numchan/2) = conv(raw_noise,HRIRNO(:,2*cont),'valid');
end  

clear cont;

%**************************************************************************
% VAD - Voice Activity Detection
%**************************************************************************

vad_signal = vadsohn(speech_hri(:,1),FAM,'a'); %a -default
if(plots==1) 
    plot(speech_hri(:,1)); 
    hold on; 
    plot(2.5e-3.*vad_signal); %Putting in the same scale (visualization)
end

%**************************************************************************
% Input SNR Generation
%**************************************************************************

% search for the initial nonzero vad sample -------------------------------
flag  = 1;    
cont1 = 0;   
ini   = 0;
while(flag)
    cont1 = cont1 + 1;
    if (vad_signal(cont1)>0)
        ini  = cont1;
        flag = 0;
    end
end
clear flag cont1;

% search for the final nozero vad sample ----------------------------------
flag  = 1;       
cont1 = tamsign+1;  
fim   = 0;
while(flag)
    cont1 = cont1 - 1;
    if (vad_signal(cont1)>0)
        fim  = cont1;
        flag = 0;
    end
end
clear flag cont1;

% epoch for SNR calculation -----------------------------------------------
vad          = zeros(tamsign,1);
vad(ini:fim) = 1;
norma        = sum(vad);

if( SNS == 0 )
    pot_speech = sum((speech_hri(:,1) .* vad).^2) / norma;
    pot_noise  = sum((noise_hri(:,1)  .* vad).^2) / norma;
else
    pot_speech = sum((speech_hri(:,numchan/2+1) .* vad).^2) / norma;
    pot_noise  = sum((noise_hri(:,numchan/2+1)  .* vad).^2) / norma;
end

 amp_noise = sqrt(pot_speech/pot_noise*(10^(-input_SNR/10)));
 noise_hri = amp_noise * noise_hri;

snr_front_left  = 10*log10( sum((speech_hri(:,1) .* vad).^2) / ...
    sum((noise_hri(:,1) .* vad).^2) );

snr_front_right = 10*log10( sum((speech_hri(:,numchan/2+1).* vad).^2) / ...
    sum((noise_hri(:,numchan/2+1) .* vad).^2) );

sum(sum(speech_hri));
sum(sum(noise_hri));

%**************************************************************************
% Signals at the ears
%**************************************************************************
input_signal = speech_hri + noise_hri;

norma = max( max( abs( [ input_signal noise_hri speech_hri ] ) ) );

input_signal = (FCT/norma) * input_signal;
speech       = (FCT/norma) * speech_hri;
noise        = (FCT/norma) * noise_hri;
clear norma;
if(plots==1) 

figure(); 
plot(speech(:,1));
hold on;
plot(noise(:,1)); title('Speech and noise - Left Hearing Aid - Front Mic');
legend('speech','noise');
figure(); 
plot(speech(:,(numchan/2+1)));
hold on;
plot(noise(:,(numchan/2+1))); title('Speech and noise - Right Hearing Aid - Front Mic');
legend('speech','noise');
end

%**************************************************************************
% Output files generation
%**************************************************************************

% Audio Files
% Speech audio from front microfones (L+R) 
aux = 0.9/max(abs([ speech(:,1) ; speech(:,(numchan/2+1))]));
audiowrite( FnamSpOut, [ speech(:,1) speech(:,(numchan/2+1))]*aux,FAM);

% Noise audio from front microfones (L+R)
aux = 0.9/max(abs([ noise(:,1) ; noise(:,(numchan/2+1))]));
audiowrite( FnamNoOut, [ noise(:,1) noise(:,(numchan/2+1))]*aux,FAM);

% Input signal audio from front microfones (L+R)
aux = 0.9/max(abs([ input_signal(:,1) ; input_signal(:,(numchan/2+1))]));
audiowrite( FnamInFront, [ input_signal(:,1) input_signal(:,(numchan/2+1))]*aux,FAM);

a=horzcat(input_signal(:,1),input_signal(:,2),input_signal(:,4),input_signal(:,5));

% Saving the Testbench file
fid = fopen( strcat(DIROUT,'TestBench.txt'), 'wt' );
switch M
    case 2 
            fprintf(fid,'%.15f\t%.15f\n',horzcat(input_signal(:,1),input_signal(:,4))');
    case 4 
            fprintf(fid,'%.15f\t%.15f\t%.15f\t%.15f\n',horzcat(input_signal(:,1),input_signal(:,2),input_signal(:,4),input_signal(:,5))');
    case 6
            %fprintf(fid,'%f\t %f\t %f\t %f\t %f\t %f\n',[input_signal(:,1),input_signal(:,2),input_signal(:,3),input_signal(:,4),input_signal(:,5),input_signal(:,6)]);
            fprintf(fid,'%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n',input_signal');
end
fclose(fid);

