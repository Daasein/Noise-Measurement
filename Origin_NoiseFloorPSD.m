%%This m file computes the autospectrum, prob density, and transfer funtions
%from a time domain signal acquired by the National Instruments DAQ. It is written using
%Matlab Version R2014a
%R. N. Miles

% 4138 calibration sqrt(1.365) volts/pascal with B&K 5935L set at gain of
% 50. This is 1.1683 volts/pascal with gain of 50. .01323 (volts/pasc)^2
%sqrt(.01323)*10^(-50/20) %volts/pasc=   0.00036373 volts/pascal 0.00035085 0.00035063 (for 4138)
%for 4138 S/N 2603134 0.00051287 volts/pascal using sigma. 0.00052092

%for 4130: 0.0011284 v/pa  0.001907  0.001915 0.0004342 0.00043331
%0.0019121 0.0019086  0.00042612  Based on sigma: 0.0023704 this is using
%the better calibrator which is at 1 kHz.  The other calibrator is at about
%980 Hz. 04/08/15 .00209 volts/pascal
% 0.0020729

addpath '\\lightning.binghamton.edu\jzhouresearch\ZhouShared\measurements\NI DAQ programs - noise measurement'

% 040815:  Hoy lab microhone % 0.0085782
%  0.0085799
clear all
close all
      [filename, pathname] = uiputfile('*.*', 'Enter a file to write analysis data');
      [filepath,fname,ext] = fileparts(fullfile(pathname,filename));%This handy command gets the file name without the extension
    if filename == 0, return, end
    filenamemat=fullfile(pathname,[fname,'.mat']); 
   cd(pathname); 

%% Define the calibration factors for each channel used.  These convert the acquired voltage to the proper units of the measurement
% And- store the gain that is applied to each channel before the signal is
% input to the DAC.  
% If a channel is unused just set it to 1.

 % calibration=[ 0.00040537  .00288 1 1 1]; %ch2 is volts/(mm/s) %The calibation factors in volts/unit.  i.e. if it is a microphone, these should be volts/pascal. Measured 032715 for the 4138-A-015
 % gain=[10^(40/20) 10^(30/20) 100 1  1 1]; %The gain applied to each channel outside the DAQ.  Note if the gain is x dB, then this should be set to 10^(x/20).
 % calibration=[ 0.000469472  1/5 1 0.00040537 1 1]; %ch2 is volts/(mm/s) %The calibation factors in volts/unit.  i.e. if it is a microphone, these should be volts/pascal. Measured 032715 for the 4138-A-015
 % gain=[10^(50/20) 1 5 10^(50/20)  1 1]; %The gain applied to each channel outside the DAQ.  Note if the gain is x dB, then this should be set to 10^(x/20).
 calibration=[ 1  1 1 1 1 1]; %ch2 is volts/(mm/s) %The calibation factors in volts/unit.  i.e. if it is a microphone, these should be volts/pascal. Measured 032715 for the 4138-A-015
 gain=[1 1 1 1  1 1]; %The gain applied to each channel outside the DAQ.  Note if the gain is x dB, then this should be set to 10^(x/20).
numchannels=2; %The actual number of channels being used
%% Define the units for each channel. This is essential to ensure each plot is labeled correctly
units={};
units{1}='volts';
units{2}='volts';
units{3}='volts';
units{4}='volts';
%%
numrecord=100; %The number of records to acquire.
numrecordwait=0;%Wait for numrecordwait records before averaging the data.
% num=16384; %num is the total number of time points per record. This should be a power of 2 for the FFT algorithm.
num=32768;
numtotal=numrecord*num;
fs=96000;
% fs = 204800;
Trecord=num/fs;%The amount of time in each record
totaltime=numrecord*num/fs; %total time in seconds.  
delfreq=1/Trecord; %The frequency spacing in Hertz.

%% This section creates legends for the plots to label each curve with a channel number
channellegend= {};
    for ichannel=1:numchannels
        channellegend{ichannel}=[fname, ' channel ',num2str(ichannel),' ',units{ichannel}];
    end
    %%
    B500=1; %Set B500=1 if using the DAC in B500.  Set it to anything else if using B224
NIAcquireSetUpDACForCalibrationScript040925
[audioarray,timeall]=s.startForeground;
maxvoltagelevel=max(abs(audioarray))

%Apply the calibration and gain settings:
for ichannel=1:numchannels
audioarray(:,ichannel)=audioarray(:,ichannel)/(calibration(ichannel)*gain(ichannel));
end
%%
for ichannel=1:numchannels
MeanSignal(ichannel)=mean(audioarray(:,ichannel));
RMSsignalTimeDomain(ichannel)=sqrt(mean((audioarray(:,ichannel)-MeanSignal(ichannel)).^2));
end
%%
delt=1/fs; 
numby2=num/2;
freq=linspace(0,fs/2,numby2);
specall=zeros(length(freq),numchannels);
ex2timeall=zeros(numchannels,1);
%%
for irec=1:numrecord
    trec=timeall(1+num*(irec-1):num+num*(irec-1));
    Trec=timeall(num+num*(irec-1))-timeall(1+num*(irec-1));
weighting1=(1-cos(2*pi*(trec-trec(1))/Trec)); %Create a Hanning window that has unit RMS for each record
weightingRMS=sqrt(mean((weighting1).^2));
weighting1=weighting1/weightingRMS;%Normalize the Hanning window so that has unit RMS for each record
%%
weighting=weighting1;
for ichannel=2:numchannels
    weighting=[weighting1 weighting];
end
%%
rtimeall=audioarray(1+num*(irec-1):num+num*(irec-1),:).*weighting;
if irec==1
    %%
    figure('name','acquired raw data')
plot(timeall,audioarray)
xlabel('time (seconds)')
ylabel('raw data')
    title(['raw time domain data rtimeall for record ',num2str(irec)])
legend(channellegend{:})
grid on
    figure('name','acquired raw data weighted')
    plot(trec,rtimeall)
    xlabel('time (seconds)')
    ylabel('rtimeall')
    title(['Hanning weighted time domain data rtimeall for record ',num2str(irec)])
    legend(channellegend{:})
    grid on
end
%%
averecord=mean(rtimeall,1);
fftalltemp=fft(rtimeall,num);  %
fftall=fftalltemp(1:numby2,:);%The frequency domain data has half the size of the time domain data.
specall=specall+abs(fftall).^2*2*delt/num; %calculate the autospectrum of the random signal. 
%                                          This will be a single-sided spectrum in units^2/Hz
%% The following calculates sigma for each record
for ichannel=1:numchannels
 ex2timeall(ichannel)=ex2timeall(ichannel)+sum(rtimeall(:,ichannel).^2,1)/num; %Mean square of each time record
sigmarecord(irec,ichannel)=sqrt(mean((rtimeall(:,ichannel)-averecord(ichannel)).^2));%This is the variance about the mean
end
TimeMeansquare(irec)=timeall(1+num*(irec-1));
%Calfactorrecord(irec,:)=sigmarecord(irec,:).*10.^(-gain/20); %The cal
%factor for each record
end
%%
numrecordused=numrecord-numrecordwait;
specall=specall/numrecordused;
rootpsd=sqrt(specall);
ex2timeall=ex2timeall/numrecordused; %disp([num2str(freq(ifreq)),'Hz OK.  Running....']);
%%
disp('mean square from time domain:')
ex2timeall
ex2specall=((specall(1,:)/2+sum(specall(2:numby2,:),1))*delfreq)'; %compute the mean square using the power spectrum by integrating over all frequencies. 
disp('mean square from frequency domain:')
ex2specall
%%
% RelativeErrorinMeanSquare=(ex2timeall-ex2specall)./ex2specall
RMSfromSpectrum=sqrt(ex2specall);
%estimate the probability density of the signal.
calcprob=1;
if calcprob==1
numbin=100;
bin(1:numbin+1,numchannels)=0;
avesignals=mean(audioarray,1);

for ichannel=1:numchannels
sigma(ichannel)=sqrt(mean((audioarray(:,ichannel)-avesignals(ichannel)).^2));%This is the variance about the mean

xmax=max(audioarray(:,ichannel));
xmin=min(audioarray(:,ichannel));
delx=(xmax-xmin)/numbin;
xbin(:,ichannel)=linspace(xmin,xmax,numbin+1);
    for i=1:num*numrecord
   ibin=round((audioarray(i,ichannel)-xmin)/delx)+1;
   bin(ibin,ichannel)=bin(ibin,ichannel)+1;
    end
    xprobdens(:,ichannel)=(xbin(:,ichannel)-avesignals(ichannel))/sigma(ichannel);
bin(:,ichannel)=bin(:,ichannel)./(delx*num*numrecord);
probdens(:,ichannel)=bin(:,ichannel)*sigma(ichannel);
sumofprobability(ichannel)=sum(bin(:,ichannel),1)*delx;
end
xgaussian=linspace(-5,5,400);
gaussian=exp(-xgaussian.^2/2)/sqrt(2*pi);%The Gaussian density.
%%
figure('name','probability density')
plot(xprobdens,probdens)
xlim([-5 5])%Plot the densities over 5 standard deviations
ylabel('prob dens')
xlabel('normalized x')
hold on
plot(xgaussian,gaussian)
legend(channellegend{:})
hold off
grid on
end
%
h=figure('name','root psd all');
loglog(freq,rootpsd)
ylabel('root psd (units/sqrt(Hz))')
xlabel('frequency (Hz)')
legend(channellegend{1:numchannels})
grid on
saveas(h,[fname,' root psd noise']);
%%
save(filenamemat, 'freq', 'specall')