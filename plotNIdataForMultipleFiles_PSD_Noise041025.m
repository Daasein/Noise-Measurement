%%
%close all
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
freq=allthedata(1).freq; %The frequency range should be the same for all of the data. Use just the first half of the data.
numfreq=length(freq);
numchannels=length(allthedata(1).specall(1,:));
numfiles=ii;
% 
for ichannel=1:numchannels; %Plot all six channels
for ifile=1:numfiles
psdplot(:,ifile)=allthedata(ifile).specall(:,ichannel);
end

% figure('name',['psd for channel ',num2str(ichannel)])
% loglog(freq,psdplot)
% xlabel('frequency (Hz)')
% ylabel('PSD (v^2/Hz)')
% title(['psd for channel ',num2str(ichannel)])
% grid on
% legend(datafiles)

figure('name',['root psd for channel ',num2str(ichannel)])
loglog(freq,sqrt(psdplot))
xlabel('frequency (Hz)')
ylabel('PSD (unit/root(Hz))')
title(['psd for channel ',num2str(ichannel)])
grid on
legend(datafiles)


end

