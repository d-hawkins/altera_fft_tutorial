% -----------------------------------------------------------------------------
% altera_burst_fft_gaussian_analysis.m
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
%
% For 18-bit data and twiddles:
%
%           |         Altera Burst           |     Radix-4 BFP
%           | ------------------------------ | --------------------
%   LF (dB) | QN (dB)  NPR (dB)   Exp Codes | QN (dB)  NPR (dB)   Exp Codes
%   ------- | -------  --------   --------- | -------  --------   ---------
%     -15   |  -88.5     73.5        7,8    |  -88.1     73.1        6,7
%     -20   |  -90.6     70.6        6,7    |  -92.1     72.1        5,6
%     -25   |  -92.1     67.1        5,6    |  -95.2     70.2        4,5
%     -30   |  -93.6     63.6        5      | -100.8     70.8        3,4
%     -35   |  -93.9     58.9        4,5    | -105.7     70.7        3
%     -40   |  -94.0     54.0        3,4    | -107.9     67.9        2,3
%     -45   |  -94.1     49.1        2,3    | -108.4     63.4        1,2
%     -50   |  -94.1     44.1        2      | -108.5     58.5        0,1
%
% For documentation figures use -15, -25, and -35dB.
%
LF_dB = -15;

% Use Altera FFT or Radix-4 BFP model?
%use_altera = 0;
use_altera = 1;

% Figure generation
if (use_altera)
	pdfbase = sprintf('altera_burst_gaussian_%d_%d', Bx, Bw);
else
	pdfbase = sprintf('radix4_bfp_gaussian_%d_%d', Bx, Bw);
end
pdfbase = '';

% -----------------------------------------------------------------------------
% FFT function handles
% -----------------------------------------------------------------------------
%
if (use_altera)
	fft_fixed = @(x) altera_burst_fft_model(x,Bx,Bw,N,0);
else
	fft_fixed = @(x) fft_radix4_bfp_model(x,Bx,Bw);
end

% MATLAB FFT
fft_double = @(x) fft(x);
%fft_double = @(x) fft_radix4_model(x,Bw);

% -----------------------------------------------------------------------------
% Analysis
% -----------------------------------------------------------------------------
%
fft_analysis_gaussian(fft_fixed, fft_double, Bx, N, M, LF_dB, pdfbase)
