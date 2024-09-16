% -----------------------------------------------------------------------------
% altera_burst_fft_model.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Altera burst FFT model.
%
% -----------------------------------------------------------------------------
% Altera Burst I/O FFT Mex
% ------------------------
%
% This script uses the mex file Sfftmodel.mexw64.
%
% The mex file can be obtained from within the Quartus installation, eg.,
%
% intel/prime_lite/22.1std/ip/altera/dsp/altera_fft_ii/src/matlab/lib
%
% Alternatively, during the Altera FFT generation process, generate the
% example design, and the mex file is located in the example directory, eg.,
%
% altera_burst_fft_example/Matlab_model
%
% Copy the mex file to the same directory as this script, or use
% MATLAB addpath.
%
% -----------------------------------------------------------------------------

function X = altera_burst_fft_model(x,Bx,Bw,N,inverse)

% Parameterization Space
THROUGHPUT=4;
ARCH=2;

% Check for Altera model
if (exist('Sfftmodel') ~= 3)
	error('Error: Please setup the Altera FFT mex file!')
end

% Altera model
[roc,ioc,eoc] = Sfftmodel(real(x),imag(x),N,THROUGHPUT,ARCH,Bx,Bw,inverse);

% Result structure
X.data     = roc + 1j*ioc;
X.exponent = -unique(eoc);
