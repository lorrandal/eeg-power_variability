clc
clear all
close all

%% LOAD DATA 
% SIGNLAS: F4C4, DELTA2, DELTA2_PV, BFV2.
load data/dataset09.mat

% BFV2 - blood flow rate (cm/sec) measured with ecodoppler at the right 
% cerebral artery with sampling frequency Fs BFV = 0.5 Hz.

% F4C4 - EEG (mV) measured from a frontal derivation located on the right 
% side of the scalp with sampling frequency Fs EEG = 512 Hz.

%%  PART 1
% Extract the delta component (DELTA2: 0-4 Hz) from the F4C4 signal using an
% elliptical low pass filter, with maximum deviation in the pass and stop 
% bands from the ideal values of 0.05 dB and a transition band of 0.5 Hz. 
% Implement the filter to eliminate phase distortion.

% Filter parameters
Rp_db = -20*log10(0.95); % passband ripple
Rs_db = -20*log10(0.05); % stopband ripple
Wp = 4/(Fs_EEG/2); % passband edge frequency      
Ws = 4.5/(Fs_EEG/2); % stopband edge frequency

[order, Wn] = ellipord(Wp, Ws, Rp_db, Rs_db);

% Filter coefficients
[b,a] = ellip(order, Rp_db, Rs_db, Wn);

% Extracting DELTA2 component
DELTA2 = filtfilt(b, a, F4C4);

%% PART 2
% Report the difference equation of the filter, plot on a 
% frequency scale Gain and Phase of its frequency response.

[H,F] = freqz(b, a, 512, Fs_EEG);

figure
subplot(211)
plot(F,abs(H))
title(['Frequency response of the Lowpass filter'])
ylabel('Gain')
xlabel('Frequency (Hz)')
xlim([0 20])
ylim([0 1.2])
subplot(212)
plot(F, angle(H))
ylabel('Phase (rad)')
xlabel('Frequency (Hz)')
xlim([0 20])

% Round the coefficients
b = roundn(b, -4);  
a = roundn(a, -4); 
disp(['Difference equation of a ' num2str(order) 'th order filter ' ])
disp('')
disp(['y(n) = ' num2str(-a(2)) '*y(n-1)' num2str(-a(3)) '*y(n-2) + ' num2str(-a(4)) '*y(n-3)' ...
    num2str(-a(5)) '*y(n-4) + ' num2str(-a(6)) '*y(n-5) + ' num2str(b(1)) '*x(n)' num2str(b(2)) '*x(n-1) + ' ...
     num2str(b(3)) '*x(n-2) + ' num2str(b(4)) '*x(n-3)' num2str(b(5)) '*x(n-4) + ' num2str(b(6)) '*x(n-5)'])

%% PART 3
% Plot on two graphs 5 seconds of signal F4C4 and 5 seconds of signal DELTA2
% indicating on the graphs with dotted lines the mean ± SD.

t = 0:1/Fs_EEG:5-1/Fs_EEG;
% Mean and std in 5 sec interval
mean_F4C4 = mean(F4C4(1:length(t)));
SD_F4C4 = std(F4C4(1:length(t)));

mean_DELTA2 = mean(DELTA2(1:length(t)));
SD_DELTA2 = std(DELTA2(1:length(t)));

% 5 sec plot with mean  and +- std
figure
subplot(211)
plot(t, F4C4(1:length(t)), t, mean_F4C4*ones(1,length(t)), '--r', t, ...
    (mean_F4C4+SD_F4C4)*ones(1,length(t)), '--r', t, ...
    (mean_F4C4-SD_F4C4)*ones(1,length(t)), '--r')
title('Signal F4C4')
ylabel('mV')
xlabel('time (s)')
ylim([-15 25])

subplot(212)
plot(t, DELTA2(1:length(t)), t, mean_DELTA2*ones(1,length(t)), '--r',t, ...
    (mean_DELTA2+SD_DELTA2)*ones(1,length(t)), '--r', t, ...
    (mean_DELTA2-SD_DELTA2)*ones(1,length(t)),'--r')
title('Signal DELTA2')
xlabel('time (s)')
ylabel('mV')
ylim([-8 5])

%% PART 4
% Segment the DELTA2 signal in 2sec intervals.

% Create a matrix with 1024 columns (512x2 -> no. samples in 2 sec)
n_col = Fs_EEG*2; % 1024 = length of the segments
n_row = length(DELTA2)/n_col; % 150 = no. segments

matr_segm = zeros(n_row, n_col);
k = 1;
for i = 1:n_row
    matr_segm(i,:) = DELTA2(k:k+n_col-1); % n_col-1023
    k = k+n_col; % n_col = 1024
end

%% PART 5 
% For each interval of 2 seconds, evaluate the spectrum using the periodogram
% method and integrate (numerically) the spectrum to evaluate the total 
% power associated with the delta band.

spectrum_segm = zeros(n_row,n_col); % matrix that contains the spectrums of the segments

% Periodogram on each segment (row)
for i = 1:n_row
    P_period = (1/n_col)*abs(fft(matr_segm(i,:)-mean(matr_segm(i,:)),n_col)).^2;
    spectrum_segm(i,:) = P_period;
end

% Take half of the samples
spectrum_segm = spectrum_segm(:,1:n_col/2);

% Frequency for plotting, take half of the samples
freq = 0:Fs_EEG/n_col:Fs_EEG-Fs_EEG/n_col;
freq = freq(1:n_col/2);

% Numerical integration performed in part 6

%% PART 6
% Construct the DELTA2 PV (Power Variability) signal which reports, at 2-second intervals (hence
% with a sampling frequency of 0.5 Hz), the total power associated with the delta band.

% Define the samplig rate for trapz
s_rate = Fs_EEG/1024;

DELTA2_PV = zeros(1, n_row);
for i = 1:n_row
    DELTA2_PV(i) = s_rate*trapz(spectrum_segm(i, 1:9)); %Integrate each segment with "trapz" and join them
end
Fs_D2PV = 0.5; %sampling frequency of DELTA2_PV


t_d2pv = 2:2:150*2; 
figure
plot(t_d2pv, DELTA2_PV, '-')
title('DELTA2 PV')
ylabel('mV^2')
xlabel('time (sec)')

%% PART 7
% Evaluate the spectrum of the DELTA2 PV signal with the direct FT method (Periodogram).

L = length(DELTA2_PV);
spectrum_D2PV = (1/L)*abs(fft(DELTA2_PV-mean(DELTA2_PV), L)).^2; % subtract the mean before computing fft

spectrum_D2PV = spectrum_D2PV(:, 1:L/2);
f_D2PV = 0:Fs_D2PV/L:Fs_D2PV-Fs_D2PV/L;
f_D2PV = f_D2PV(1:L/2);

figure
plot(f_D2PV, spectrum_D2PV)
title('Spectrum of DELTA2 PV (Periodogram) ')
ylabel('mV^2/Hz')
xlabel('Frequency (Hz)')

%% PART 8
% Identify the optimal order AR model on the DELTA2 PV signal 
% (search for the order in the range1-20)

orders = 1:20;
no = length(orders);
FPE = zeros(no, 1);
AIC = zeros(no, 1);
Jmin = zeros(no, 1);

% Calculate for each model order: MSE, Akaike's Final Prediction Error,
% Akaike Information Criterion
for ord = 1:no
    th = ar(DELTA2_PV-mean(DELTA2_PV), ord, 'yw');
    Jmin(ord) = th.noisevariance;
    FPE(ord) = fpe(th); 
    AIC(ord) = aic(th); 
end

figure
plot(orders, Jmin, 'b-o')
legend('MSE')
ylabel('MSE')
xlabel('Order')

figure
plot(orders, AIC, 'r-^')
legend('AIC')
ylabel('AIC')
xlabel('Order')

figure
plot(orders, FPE, 'g-*')
legend('FPE')
ylabel('FPE')
xlabel('Order')

ord_aic = find(AIC==min(AIC));
ord_fpe = find(FPE==min(FPE));

% Evaluation of Jmin ratio for subsequent orders
epsilon = 0.001;
p = 0;
ratio = 0;
while((ratio<=(1-epsilon))&(p<no))
    p = p+1;
    ratio(p) = Jmin(p+1)/Jmin(p);
end
ord_j = p;

disp(['Optimal order by AIC: ' num2str(ord_aic)])
disp(['Optimal order by FPE: ' num2str(ord_fpe)])
disp(['Optimal order by J: ' num2str(ord_j)])

% ANDERSON TEST
ord_test = 1;
s = DELTA2_PV-mean(DELTA2_PV); %subtract the mean from DELTA2_PV
test = ar(s, ord_test, 'yw');
err = filter(test.a, 1, s);
% err=err-mean(err);

% calculate ro(k)=R(k)/R(0)
l_err = length(err);
k_max = 50;
ro = zeros(k_max+1 ,1);
for k = 0:k_max
    temp = 0;
    for n = 0:(l_err-k-1)
        temp = temp+err(n+1)*err(n+1+k);
    end
    ro(k+1) = temp;
end
ro = ro/ro(1);
ro(1) = []; % ro(1)=1, we are interested in the other values
std_err = sqrt(1/l_err);
alpha = 0.05; % set the threshold
bb = norminv([alpha/2 1-alpha/2], 0, std_err);
beta = bb(2);
exceed = length(find(abs(ro)>beta)); % Search for values that exceed the threshold
disp(['Values that exceed the threshold are ', num2str(exceed/length(ro)*100),'%']);

% PLOT
figure
box on, hold on
plot(ro,'o--')
plot([0 length(ro)], [0 0],'r--')
plot([0 length(ro)], [bb(1) bb(1)], 'k--')
plot([0 length(ro)], [bb(2) bb(2)], 'k--')
ylim([-0.2 0.2])
title('Estimation of \rho coefficients for the Anderson test')
xlabel('k'), ylabel('Normalized values')

%% PART 9
% Difference equation of the optimal order AR model.

disp(['Difference equation of the AR model of order ' num2str(ord_test)])
disp(['y(n)= ' num2str(roundn(test.a(2),-4)) '*y(n-1) + x(n)'])

%% PART 10
% Use the optimal order AR model to estimate the spectrum of DELTA2 PV signal.

model = ar(s, 1, 'yw');
[H,f] = freqz(1, model.a, L, Fs_D2PV);
PSD = ((abs(H)).^2)*model.noisevariance;
figure
plot(f, PSD)
title('Power Spectral Density with AR model');
ylabel('mV^2/Hz')
xlabel('Frequency (Hz)')

%% PART 11
% Evaluate the Coherence function between DELTA2 PV and BFV2 signals.

% [Coer, F1] = mscohere(DELTA2_PV, BFV2, L/15 ,L/30 , L, Fs_BFV); % Trying
% several windows
% [Coer, F1] = mscohere(DELTA2_PV, BFV2, L/3, 0.5*(L/3), L, Fs_BFV);
[Coer, F1] = mscohere(DELTA2_PV, BFV2, L/5, 0.5*(L/5), L, Fs_BFV);
figure
plot(F1, Coer)
title('Coherence function between DELTA2 PV and BFV2')
xlabel('Frequency (Hz)')

%% PART 12
% Spectrum of DELTA2 PV, Periodogram vs AR model

figure
plot(f_D2PV, spectrum_D2PV, f, PSD, 'r')
legend('Periodogram','AR')
title('Spectrum of DELTA2 PV: Periodogram vs AR')
ylabel('mV^2/Hz')
xlabel('Frequency (Hz)')














