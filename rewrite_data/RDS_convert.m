%% This script reads in a Bscope frame (*.bfr) and plots it 
%  
%  - the image data is contained in binary format
%
%
clear all; close all; %clc;
fprintf('Select the first Bscope file (*.brf)');  

%% Open starting frame (*.bfr) and read according to written format (ABCDEFGH) 
[Filename, Pathname] = uigetfile('*.bfr');  	% get *.bfr file
Filebasei = [Filename(1:19)];                   % File basename with indicator
addpath(Pathname);                              % Include path for reading mfile
addpath(Pathname);                              %
eval(['Prm_' Filebasei]);                       % Open associated parameter file
% first frame #
n1 = Filename(21:24);
n1 = str2num(n1);

%% total number of frames to open
tn = inputdlg('Enter number of frames to open:');
tn = str2num(tn{1});

Ang_endpt = 0;

for n = n1+10000:tn-1+n1+10000
    s = sprintf('loading frame #%d', n);
    disp(s);
 
    % open
    Filename = sprintf('%s%d.bfr', Filebasei, n);
    Filename(length(Filename)-8)='_';
    FilePointF = [Pathname, Filename];
    fdF = fopen(FilePointF, 'rb');
    
    % -------------------------------------------------------------------------
    % read
    Az_bins= fread(fdF, 1, 'int16');         	% A read    # Sweeps in the current Image Frame
    Az_cnt = fread(fdF, Az_bins, 'int16');  	% B read    Azimuth count value
    RP_dBm = zeros(fft_bins, Az_bins);      	% Initialize    
    for mm = 1:Az_bins
        RP_dBm(:,mm) = fread(fdF, fft_bins, 'single');	% C Image data
    end
    
    fclose(fdF);
    
    Ang_plot = Az_cnt .* 0.02;
    if (Ang_endpt < max(abs(Ang_plot)))
        Ang_endpt = max(abs(Ang_plot));
    end
    
    frameWriteName = [Pathname,Filename];
    frameWriteName(end-8) = '';
    frameWriteName(end-2:end) = 'dat';
    uwEachWrite(frameWriteName, Ang_plot, RP_dBm);
end

baseWriteName = [Pathname, 'Prm_', Filebasei, '.dat'];
uwBaseWrite(baseWriteName, Ang_endpt, MagMax, MagMin, fft_dR);