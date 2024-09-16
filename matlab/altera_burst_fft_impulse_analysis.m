% -----------------------------------------------------------------------------
% altera_burst_fft_impulse_analysis.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Altera burst I/O FFT analysis.
%
% The Altera FFT block floating-point implementation causes the impulses
% to have an exponent value of 5 or 6, i.e., 5 or 6 bits of the input are
% removed in the output. The input bit-width should be 16-bits or more to
% ensure at least 10-bits are used in the output.
%
% There is no bit-growth for an impulse through an FFT, so the expected
% exponent from the Altera FFT was 3, i.e., the 3-bits that the exponent
% shifts to ensure no overflow through a radix-4 engine pass. An exponent
% code of 5 or 6 was not expected.
%
% See altera_burst_fft_constant_analysis.m for additional analysis of
% the unexpected exponent scaling.
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

% Use Altera FFT or Radix-4 BFP model?
%use_altera = 0;
use_altera = 1;

% Figure generation
if (use_altera)
	pdfbase = sprintf('altera_burst_impulse_%d_%d', Bx, Bw);
else
	pdfbase = sprintf('radix4_bfp_impulse_%d_%d', Bx, Bw);
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
fft_analysis_impulse(fft_fixed, fft_double, Bx, N, pdfbase)
