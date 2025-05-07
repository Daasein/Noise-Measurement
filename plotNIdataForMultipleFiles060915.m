%%
close all
clear all
%This M file reads data created by NIMatlabspectrumanalyzer032215.m

%Pick the files and load the data
ii = 0;
while 1
    [filename, pathname] = uigetfile('*.mat', 'pick a file to read. Click cancel to end.','MultiSelect','on');
    if ischar(filename)
        filename = {filename};
    elseif isnumeric(filename) & filename==0
        break
    end
    for iii=1:length(filename)
        ii = ii + 1;
        datafiles{ii} = filename{iii};
        pathnames{ii} = pathname;
        allthedata(ii)=load(fullfile(pathnames{ii},datafiles{ii}));
        cd(pathname);
    end
end

%% Plot data from each file
freq=allthedata(1).alldata.freq; %The frequency range should be the same for all of the data. Use just the first half of the data.
numfreq=length(freq);
numchannels=length(allthedata(1).alldata.psd(1,:));
numfiles=ii;
% 
for ichannel=1:numchannels; %Plot all six channels
for ifile=1:numfiles
psdplot(:,ifile)=allthedata(ifile).alldata.psd(:,ichannel);
transferfnplot(:,ifile)=allthedata(ifile).alldata.tranfn(:,ichannel)./allthedata(ifile).alldata.tranfn(:,1);
sensitivityplot(:,ifile)=20*log10(abs(allthedata(ifile).alldata.tranfn(:,ichannel)./allthedata(ifile).alldata.tranfn(:,1))/10);%Sensitivity for Knowles
coherenceplot(:,ifile)=allthedata(ifile).alldata.coherence(:,ichannel);
% OASPLplot(:,ifile)=allthedata(ifile).alldata.OASPLrec(:,ichannel);
% timeOASPL(:,ifile)=allthedata(ifile).alldata.timerec;
end

figure('name',['psd for channel ',num2str(ichannel)])
loglog(freq,psdplot)
xlabel('frequency (Hz)')
ylabel('PSD (v^2/Hz)')
title(['psd for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['root psd for channel ',num2str(ichannel)])
loglog(freq,sqrt(psdplot))
xlabel('frequency (Hz)')
ylabel('PSD (v/root(Hz))')
title(['psd for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['coherence for channel ',num2str(ichannel)])
semilogx(freq,coherenceplot)
xlabel('frequency (Hz)')
ylabel('coherence')
title(['coherence for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['linear freq transfer function magnitude Re ch 1 for channel ',num2str(ichannel)])
semilogy(freq,abs(transferfnplot))
xlabel('frequency (Hz)')
ylabel('H1 transfer function (v/v)')
title(['transfer function magnitude for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['transfer function magnitude  Re ch 1 for channel ',num2str(ichannel)])
loglog(freq,abs(transferfnplot))
xlabel('frequency (Hz)')
ylabel('H1 transfer function (v/v)')
title(['transfer function magnitude for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['transfer function phase Re ch 1 for channel ',num2str(ichannel)])
semilogx(freq,angle(transferfnplot)*180/pi)
xlabel('frequency (Hz)')
ylabel('H1 phase (degrees)')
title(['transfer function phase for channel ',num2str(ichannel)])
grid on
legend(datafiles)

figure('name',['sensitivity in dB Re 1 v per .1 pa for channel ',num2str(ichannel)])
semilogx(freq,sensitivityplot)
xlabel('frequency (Hz)')
ylabel('dB Re 1 v per .1 pa')
title(['sensitivity for channel ',num2str(ichannel)])
grid on
legend(datafiles)
% 
% figure('name',['OASPL for channel ',num2str(ichannel)])
% plot(timeOASPL,OASPLplot)
% xlabel('time (seconds)')
% ylabel('OASPL (dB)')
% title(['OASPL for channel ',num2str(ichannel)])
% grid on
% legend(datafiles)

end

