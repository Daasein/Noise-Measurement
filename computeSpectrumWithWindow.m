function [freq, rootPSD] = computeSpectrumWithWindow(signal, fs, plotFlag, averageCount)
% computeSpectrumWithWindow - Multi-channel spectrum with Hanning window and averaging
%
% Syntax:
%   [freq, rootPSD] = computeSpectrumWithWindow(signal, fs, plotFlag, averageCount)
%
% Inputs:
%   signal        - Time-domain signal [N x numChannels]
%   fs            - Sampling frequency in Hz
%   plotFlag      - (optional) true/false for plotting (default: true)
%   averageCount  - (optional) Number of segments to average (default: 1)
%
% Outputs:
%   freq          - Frequency vector [Hz]
%   rootPSD       - Root PSD [numFreqs x numChannels]

if nargin < 3, plotFlag = true; end
if nargin < 4, averageCount = 1; end

signal = double(signal);
[N, numChannels] = size(signal);

segmentLength = 2^floor(log2(N / averageCount));
numby2 = segmentLength / 2;
delt = 1 / fs;
freq = linspace(0, fs/2, numby2)';
specSum = zeros(numby2, numChannels);

for i = 1:averageCount
    idxStart = (i-1)*segmentLength + 1;
    idxEnd = idxStart + segmentLength - 1;
    if idxEnd > N
        warning('Incomplete segment skipped at end of signal.');
        break;
    end
    segment = signal(idxStart:idxEnd,:);

    % Apply Hanning window with unit RMS
    t = (0:segmentLength-1)' / fs;
    w = 1 - cos(2 * pi * (t - t(1)) / (t(end) - t(1)));
    w = w / sqrt(mean(w.^2));  % Normalize RMS
    windowed = segment .* w;

    % FFT and autospectrum
    fftSeg = fft(windowed, segmentLength);
    fftOneSide = fftSeg(1:numby2,:);
    spec = abs(fftOneSide).^2 * 2 * delt / segmentLength;

    specSum = specSum + spec;
end

specAvg = specSum / averageCount;
rootPSD = sqrt(specAvg);

% Optional plotting
if plotFlag
    figure;
    loglog(freq, rootPSD);
    xlabel('Frequency (Hz)');
    ylabel('Root PSD (unit/\surdHz)');
    title(['Averaged Spectrum (', num2str(averageCount), ' segments)']);
    grid on;
end
end
