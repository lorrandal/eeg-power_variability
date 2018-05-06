# eeg-power_variability

This project concerns the combined study of two signals: EEG power variability (EEG-PV) and middle cerebral artery blood 
flow velocity (BFV).
A first phase of the study was conducted in subjects under anaesthesia, so in a non-physiological situation where 
heartbeat and breathing oscillations are limited. The results were published in (*Zanatta et al., The human brain 
pacemaker, Neuroimage 2013*). The study is now extended to subjects studied under normal conditions.
The dataset used was provided in the Biological Signals Processing course, held by Professor Toffolo, Bioengineering,
Department of Information Engineering, University of Padua, academic year 2014/2015.

The two signals provided in the `data` folder are:

`BFV2` : blood flow rate (cm/sec) measured with ecodoppler at the right 
cerebral artery with sampling frequency Fs BFV = 0.5 Hz.

`F4C4` : EEG (mV) measured from a frontal derivation located on the right 
side of the scalp with sampling frequency Fs EEG = 512 Hz.

To extract the delta component (`DELTA2`: 0-4Hz) from `F4C4` a **low pass filter** was used, the optimum order of 5 was found with 
the `ellipord` function while the `b` and `a` parameters were obtained thanks to the `ellip` function. To eliminate phase distortion,
Forward-Backward filtering was used using the `filtfilt` function.

The **Difference equation** of the 5th order filter is:

```y(n) = 4.9382y(n-1)-9.7584y(n-2) + 9.6455y(n-3)-4.7689y(n-4) + 0.9435y(n-5) + 0.0044x(n)-0.0131x(n-1) + 0.0087x(n-2) + 0.0087x(n-3)-0.0131x(n-4) + 0.0044x(n-5)```

*Gain* and *Phase* of the **Frequency Response** of the filter are shown in Fig. 1:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot1.svg)

The result of the filtering process is `DELTA2`, shown in Fig 2.:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot2.svg)

Once the delta component (`DELTA2` signal) had been obtained through filtering, it was segmented into 150 segments 
by 1024 samples at 2-second intervals. To construct the **Power Variability** signal (`DELTA2_PV`), the spectrum of each segment
was calculated using the **Periodogram** method, having first removed the average for each segment. Then each spectrum was 
integrated, with `trapz`, to obtain the power relative to the delta band. `DELTA2_PV` was formed by combining the values for 
each segment.

The obtained signal `DELTA2_PV` is shown in Fig. 3:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot3.svg)

The spectrum of `DELTA2_PV` signal, obtained by the Periodogram method is shown in Fig 4:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot4.svg)

Subsequently the **AR model** of optimal order was identified on the `DELTA2_PV` signal.
To find the optimal, three indicators were evaluated for every order between [1:20]: **MSE**, **Akaike's Final Prediction 
Error (FPE)**, **Akaike Information Criterion (AIC)**. AIC and FPE returned 1 as  optimal order while MSE returned 3.

The results for MSE, AIC and FPE are shown in Fig. 5, Fig. 6, Fig. 7:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot5.svg)


![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot6.svg)

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot7.svg)

**Anderson–Darling test** was carried out to verify the whiteness of the prediction error and the level of significance
was set at α = 5%.

The results of Anderson test are show in Fig. 8:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot8.svg)

Accordingly, by adopting order 1, the **difference equation** of the **AR** model is:

```y(n)= -0.0664*y(n-1) + x(n)```


Then, the optimal order AR model is used to estimate the spectrum of DELTA2 PV signal.
The model cannot explain the data properly and seems to do oversmoothing.
The **Power Spectral Density** obtained in this way is shown in Fig. 9:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot9.svg)

As a consequence of what has just been said, the **Spectral Coherence** is calculated between 
`DELTA2_PV` obtained by the Periodogram's method and `BFV2`.

The `mscohere` function was used to calculate the Coherence function between `DELTA2_PV` and `BFV2`.
After a tuning phase, a window containing 30 samples (L/5) was used, where `L = 150` is the number of `DELTA2_PV`
and `BFV2` samples. As overlap a number of samples equal to 50% of the window was used. 
The Spectral Coherence has a maximum of `max = 0.5911` at a frequency of 0.095 Hz indicating the presence of a
possible casuality link between `DELTA2_PV` and `BFV2`.

Coherence function is shown in Fig. 10:

![alt text](https://github.com/lorrandal/eeg-power_variability/blob/master/plot10.svg)




