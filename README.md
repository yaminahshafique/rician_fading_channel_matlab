

% % Note: This code was implemented on the R2023a version of MATLAB and
% % uses the 8.0 version of Communications Toolbox. 

clear all

format long

%Defining Parameters of Input Signal
fc = 1100e06; 
bit_rate = 900e06;
samp_per_bit = 32;
fs = bit_rate*samp_per_bit; %Sample rate
symbol_rate = bit_rate/4; %Symbol Rate = Bit Rate/4 because each symbol has 4 bits
symbol_period = 1/(symbol_rate); %Symbol Period

sig = readmatrix("16-QAM_1100 Mhz.txt"); %Loading the input signal into a matrix

l1 = length(sig); %Length of Signal

t = sig(:,1); %Storing the Time Values
a = sig(:,2); %Storing the Amplitude Values
a_max = max(a); %Maximum Value in the signal
a_min = min(a); %Minimum Value in the signal

%Plotting the QAM Signal
figure; 
plot(t,a,'b');
xlabel("Time (seconds)"); ylabel("Amplitude"); title("Plot of 16-QAM Signal vs. Time");

% Used in report:
% % figure;
% plot(t(1:2000),a(1:2000),'b');
% xlabel("Time (seconds)"); ylabel("Amplitude"); title("Plot of 1st 2000 samples 16-QAM Signal vs. Time");


%%% Details of the Signal:
%%% The signal provided to us is Quadrature Amplitude Modulated Signal with a
%%% modulation order of 16 and carrier frequency of 1100 MHz. Modulation
%%% order of 16 implies that for each symbol there will 4 bits. This number
%%% is obtained by using the following formula
bits_per_symbol = log2(16);
fprintf("The number of bits per symbol in the signal is %d",bits_per_symbol)
%%% Further the ideal 16 QAM signal consists of 16 distinct symbols, each
%%% representing a distinct signal level. The QAM signal is also
%%% considered to be a combination of AM and PSK because information is encoded
%%% in both amplitude and phase of the QAM signal.



%Implementation of Rician Fading Channel

%1-Determining/Defining the Channel Parameters:

% % Determination of Maximum Doppler Shift Frequency % %
% Assuming a mobile velocity of 30 m/s, the maximum doppler shift is given
% by:

v = 30; %Mobile Speed
c = 3e08; %Speed of Light in m/s
maxDopplerShift = (v*fc)/c
K = 10; %Based on MATLAB website example
pathGains = [0,-3,-6,-9]; %in dB
pathDelays = [0,0.2e-06,0.4e-06,0.8e-06]; %Path Delays for Frequency-selective fading

n = length(pathGains); %Length of Path Gains = Length of Path Delays

%The path gains are in fact the losses associated with each path.
% As the delay increases, so do the losses encountered on a certain path.
% It is, for this this reason, that the path gain is increasing negatively with delays.
% This can be visualized in the following plot:

figure;
plot(pathDelays,pathGains,'b-*'); xlabel('Path Delays (seconds)'); ylabel('Path Gains (dB)'); title('Variation in Path Gain with Path Delay')

%2-Defining the Channel

% There are two important cases to consider here:
% 1. Wideband/Frequency-Selective Fading: The received signal envelope has bumps/deep fades at certain points.
% This happens when the delay span is significantly larger than the symbol
% period of the transmitted signal and power levels of the signal vary with frequency. In such a case, intersymbol
% interference (ISI) occurs. To understand the effects of frequency
% selective fading, we have set the total delay span such that it is
% approximately 100 times larger than the symbol period.

% 2. Narrowband/Flat Fading: In this case, the received signal envelope is
% almost flat. This happens when the delay span is much smaller than the
% symbol period. There is little to no ISI. To observe this phenomenon, we
% have reduced the bit rate of the original signal and hence increased the
% symbol period. The best approach would be to regenerate the same signal with
% reduced bit keeping all other parameters same. However, with available
% resources, we have used what may be called a naive approach. This isn't
% the best way to observe the phenomenon, however provides satisfactory
% results for our analysis. For reducing the bit rate, we have simply downsampled the
% original signal by a factor of 10 and stored it in a different matrix. The downsampling factor
% of 10 was chosen so that bit rate also gets reduced by a factor 10. We
% have then passed this new signal through the channel.


%To implement the channel, we have used the built-in MATLAB command and set
%the parameters based on an example on the MathWorks website. In order to
%make sure that the NLOS random component is not random, we have set the
%'RandomStream' parameter to 'mt19937ar' which sets the seed value to 73.
%This will make sure that the simulator uses the same seed everytime the
%code is run. It is important to note here that the NLOS component is
%essentially random in practical settings and has been "fixed" only for the
%sake of analysis. Removing the 'RandomStream' argument or setting its
%value to 'Global stream' will generate a different signal everytime.


sig_downsampled = downsample(sig,10); %Downsampling the signal
a_downsampled = sig_downsampled(:,2); %New amplitude values
t_downsampled = sig_downsampled(:,1); %New time values
tsamp = t_downsampled(2,1) - t_downsampled(1,1); %Sampling period
% fs_downsampled = 1/tsamp; %New Sampling frequency
% new_bit_rate = fs_downsampled/samp_per_bit;
% new_symbol_rate = new_bit_rate/4; 
% new_symbol_period = 4/new_bit_rate;

l2 = length(a_downsampled); %Length of downsampled signal

figure;
plot(t_downsampled,a_downsampled,'b');
xlabel("Time (seconds)"); ylabel("Amplitude"); title("Plot of Down sampled QAM signal vs. Time");

% % Used in report
% % adsamp1 = a_downsampled(1:2000);
% % tdsamp1 = t_downsampled(1:2000);
% % figure;
% % plot(tdsamp1,adsamp1,'b');
% % xlabel("Time (seconds)"); ylabel("Amplitude"); title("Plot of 1st 2000 Samples of 16-QAM Signal vs. Time with Reduced Bit Rate");

ricianChannel1 = comm.RicianChannel('SampleRate',fs,'KFactor',K,'AveragePathGains',pathGains,'PathDelays',pathDelays, ...
    'MaximumDopplerShift',maxDopplerShift,'DopplerSpectrum',doppler('Flat'), ...
    'RandomStream','mt19937ar with seed');

% The same sampling frquency as the original signal has been passed because
% otherwise the channel would have changed altogether and the results
% would'nt be as per our expectations.


release(ricianChannel1);


freq_selective_faded_sig = ricianChannel1(a); %Finding the Output Signal
flat_faded_sig = ricianChannel1(a_downsampled); %Finding the Output Signal


figure;
plot(t,freq_selective_faded_sig,'b');
xlabel('Time (seconds)'); ylabel('Amplitude'); title('Rician Faded Signal (Frequency Selective Fading)');

figure;
plot(t_downsampled,flat_faded_sig,'b');
xlabel('Time (seconds)'); ylabel('Amplitude'); title('Rician Faded Signal (Flat Fading)');



%3-EVM Calculation
%To perform EVM calculation, we have used the built-in MATLAB comm.EVM
% In order to find the EVM values, we will use a loop
%to find the received signal for all the paths individually

evm = comm.EVM(MeasurementIntervalSource = 'Custom with periodic reset', Normalization='Average reference signal power');
release(evm);


%For frequency-selective fading:
rec_sig1 = zeros(l1,n); %Defining a matrix to store all the received signals
rms_evm_vals1 = zeros(1,n); %Defining a vector to store all the rms EVM values

%For flat fading:
rec_sig2 = zeros(l2,n); %Defining a matrix to store all the received signals
rms_evm_vals2 = zeros(1,n); %Defining a vector to store all the rms EVM values


%Finding the Received Signal for every path
for k = 1:n
    
ricianChannel1 = comm.RicianChannel('SampleRate',fs,'KFactor',K,'AveragePathGains',pathGains(k),'PathDelays',pathDelays(k), ...
    'MaximumDopplerShift',maxDopplerShift,'DopplerSpectrum',doppler('Flat'), ...
    'RandomStream','mt19937ar with seed','NormalizePathGains',true);

release(ricianChannel1);



rec_sig1(:,k) = ricianChannel1(a);
rec_sig2(:,k) = ricianChannel1(a_downsampled);


end

%Using the round function to find closest integer to the original values as
%well using using the abs function to find the absolute value because the
%evm function doesn't take any fractional/complex numbers in its arguments


rx1 = round(abs(rec_sig1));
rx2 = round(abs(rec_sig2));

tx1 = round(abs(a));
tx2 = round(abs(a_downsampled));

% Loop to find EVM values for all the paths
for x = 1:n
rms_evm_vals1(x) = evm(tx1,rx1(:,x));
rms_evm_vals2(x) = evm(tx2,rx2(:,x));
end




% %Calculation of Distance Values
% In order to calculate the distance values, we have used the formula 
% Distance = Speed*Time where Speed = Speed of Light and 
% Time = Time Delays for the individual paths
distance = c.*pathDelays; %Distance in meters


%4-Plotting EVM Values

figure;
plot(distance,rms_evm_vals1,'b-*'); %Frequency selective fading
xlabel('Distance (m)'); ylabel('EVM Values (%)'); title('Plot of EVM Values vs. Distance');

hold on

plot(distance,rms_evm_vals2,'r-o'); %Flat fading
xlabel('Distance (m)'); ylabel('EVM Values (%)'); title('Plot of EVM Values vs. Distance');
legend('Frequency-selective Fading', 'Flat Fading');


%Signal Power
%Ideally, the power of the frequency-selective faded signal should be
%lower than that of the flat-faded signal because the in case of flat
%fading, the entire signal is faded by approximately the same factor and there is little to no
% inter symbol interference. The same isn't true for frequency selective fading
% because different parts of the signal aren't faded by the same factor.


%Note: Power is in Watts
transmit_power1 = mean((abs(a).^2)) %Average power of original transmitted signal
transmit_power2 = mean((abs(a_downsampled).^2)) %Average power of downsampled transmitted signal
freq_selective_faded_sig_power = mean(abs(freq_selective_faded_sig).^2) %Average Power of frequency selective faded signal
flat_faded_sig_power = mean(abs(flat_faded_sig).^2) %Average Power of flat-faded signal power








