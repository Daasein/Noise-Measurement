%   NIAcquireSinusoidSetUpDAC040719
%
%fs=1/(time(2)-time(1));%The sampling rate based on the input time array
%%
%create sinusoidal signals that step through frequencies from freqmin to
%freqmax with nstep frequency steps. Create ncycles cycles for each step
% Note that fs is the sampling rate of the DAQ.
%
% This script requires the array Lowpass for the crossover
%
%%
% num=length(signalf);
% totaltime=timef(num);
%% This section is from the Matlab data acquisition toolbox
dd=daq.getDevices;
s=daq.createSession ('ni');
s.Rate=fs; %The sampling rate in Hertz
s.DurationInSeconds=totaltime; %The duration of the sample in seconds.
% s.addAnalogInputChannel('PXI1Slot3', 0,'voltage'); %Acqure data on all 6 input channels
% s.addAnalogInputChannel('PXI1Slot3', 1,'voltage'); %Acqure data on all 6 input channels
% s.addAnalogInputChannel('PXI1Slot2', 0,'voltage');
% s.addAnalogInputChannel('PXI1Slot2', 1,'voltage');
% s.addAnalogInputChannel('PXI1Slot2', 2,'voltage');
% s.addAnalogInputChannel('PXI1Slot2', 3,'voltage');
s.addAnalogInputChannel('PXI1Slot3', 0,'voltage'); %set up the desired number of channels
if numchannels > 1
    s.addAnalogInputChannel('PXI1Slot3', 1,'voltage'); 
end
if numchannels > 2
    s.addAnalogInputChannel('PXI1Slot2', 0,'voltage');
    % Modify 04/30/25, try to add another daq device
    % s.addAnalogInputChannel('Dev1', 0,'voltage');
end
if numchannels > 3
    s.addAnalogInputChannel('PXI1Slot2', 1,'voltage');
%    Try this to add a pcb accelerometer in the 3rd analog input channel
%    s.addAnalogInputChannel('PXI1Slot2', 1, 'Accelerometer'); 
%    s.Channels(1).Sensitivity = 0.00922; Add the accelerometer sensitivity
%    in volts/g here.  
end
if numchannels > 4
    s.addAnalogInputChannel('PXI1Slot2', 2,'voltage');
end
if numchannels > 5
    s.addAnalogInputChannel('PXI1Slot2', 3,'voltage');
end
% s.addAnalogOutputChannel('PXI1Slot3', 0,'voltage'); %Send data on two output channels
% s.addAnalogOutputChannel('PXI1Slot3', 1,'voltage'); %Send data on two output channels
% 
if numchannels>2
% addTriggerConnection(s,'PXI1Slot2/PXI_Trig0','PXI1Slot3/PXI_Trig0','StartTrigger');%This is essential to synchronize the channels
end
for ichannel=1:numchannels
s.Channels(ichannel).Range=[-0.01,0.01];%Set the ranges of all channels
if ichannel == 3
    continue
end
s.Channels(ichannel).Coupling='AC';
% s.Channels(ichannel).TerminalConfig
end
% numave=10;
% data0=signalf;
% noisemax=max(abs(data0));
% data0=data0*maxvoltagelevel/noisemax;%* volt maximum signal = maxvoltagelevel volts
% data1=data0;
