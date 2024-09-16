% -----------------------------------------------------------------------------
% altera_burst_fft_npr_analysis.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Altera burst I/O FFT analysis.
%
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------
%
% Number of samples
N = 1024;

% Input bit-width
Bx = 18;

% Twiddle bit-width
Bw = 18;

% Number of estimates
M = 1000;

% Gaussian signal input loading factor
LF_dB = [-70 -60 -50 -40 -35 [-30:2:-10] -5 0];

% Use Altera FFT or Radix-4 BFP model?
%use_altera = 0;
use_altera = 1;

% Figure generation
if (use_altera)
	pdfbase = sprintf('altera_burst_npr_%d_%d', Bx, Bw);
else
	pdfbase = sprintf('radix4_bfp_npr_%d_%d', Bx, Bw);
end
pdfbase = '';

% -----------------------------------------------------------------------------
% FFT function handles
% -----------------------------------------------------------------------------
%
if (use_altera)
	fft_fixed = @(x) altera_burst_fft_model(x,Bx,Bw,N,0);
	fft_fixed_name = 'Altera FFT';
else
	fft_fixed = @(x) fft_radix4_bfp_model(x,Bx,Bw);
	fft_fixed_name = 'Radix-4 BFP FFT';
end

% MATLAB FFT
fft_double = @(x) fft(x);
%fft_double = @(x) fft_radix4_model(x,Bw);

% -----------------------------------------------------------------------------
% Analysis
% -----------------------------------------------------------------------------
%
fft_analysis_npr(fft_fixed, fft_double, Bx, N, M, LF_dB, fft_fixed_name, pdfbase)