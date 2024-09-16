% -----------------------------------------------------------------------------
% altera_burst_fft_exp_codes_analysis.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Altera burst I/O FFT analysis.
%
% The Altera implementation of block floating-point has two implementation
% details that are not documented.
%
% 1. The input is unconditionally scaled by 2-bits
%
%    The Altera burst I/O FFT appears to scale the input by 2-bits regardless
%    of whether that is necessary or not. This is likely due to the Altera
%    FFT implementing the block exponent scaling after the radix-4 stage,
%    and hence the input data has to pass through the radix-4 twice before
%    the data can be scaled for the previous time through the radix-4. That
%    means the radix-4 stage needs 5 guard bits rather than 3 guard bits.
%
%    The shift is not unconditionally required. The core could check the
%    2 MSBs of each word as they are loaded into RAM to determine if scaling
%    is required after the first radix-4 pass.
%
% 2. Positive and negative numbers result in different scaling
%
%    The cause of this discrepancy appears to be due to Altera's choice to
%    round up positive values, rather than round toward zero.
%
%    a) Full-scale positive input impulse
%
%       >> x = zeros(1024,1);x(1)=2^(Bx-1)-1;
%       >> X=altera_burst_fft_model(x,Bx,Bx,N,0)
%            data: [1024×1 double]
%        exponent: 6
%
%        2^(Bx-1)-1 = 1FFFFh with an exponent code of 6, implies a
%        binary point at 01_1111_1111_11.11_1111 and the 12-bit MSBs round up
%        to the 13-bit number 0_1000_0000_0000 (800h), where an additional bit
%        was required to represent the sign. The 18-bit number using during
%        the next radix-4 pass is 00_0000_1000_0000_0000 (00800h) or 2048.
%
%    b) Negative of full-scale positive input impulse
%
%       >> x = zeros(1024,1);x(1)=-(2^(Bx-1)-1);
%       >> X=altera_burst_fft_model(x,Bx,Bx,N,0)
%            data: [1024×1 double]
%        exponent: 5
%
%        -(2^(Bx-1)-1) = 20001h with an exponent code of 5, implies a
%        binary point at 10_0000_0000_000.0_0001 and the 12-bit MSBs round
%        to the 13-bit 1_0000_0000_0000 (1000h), where the sign bit is still
%        correct. The 18-bit number using during the next radix-4 pass is
%        11_1111_0000_0000_0000 (3F00h) or -4096.
%
%     Altera's choice to round-up the positive value results in the sign bit
%     needing to use one more bit than would have been needed if the positive
%     fraction had been truncated (rounded toward zero).
%
%     For example, consider two cases where the block scaling logic detects
%     that the 3 guard bits of an 18-bit number do not match the sign bit:
%
%     i) 01_1111_1111_1111_1111
%         * Divide-by-8 to replicate the sign bit three more times:
%           01_1111_1111_1111_1.111
%         * Round-towards-zero:
%           00_0011_1111_1111_1111
%         * There are now 3 guard bits in the MSBs
%         * floor(131071/8) = 16383 (0x03FFF)
%
%     ii) 10_0000_0000_0000_0000
%         * Divide-by-8 to replicate the sign bit three more times:
%           10_0000_0000_0000_0.000
%         * Round-towards-zero:
%           11_1100_0000_0000_0000
%         * There are now 3 guard bits in the MSBs
%         * floor(-131072/8) = -16384 (0x3C000)
%
%     What is the impact of this choice? It causes a signal-to-noise loss
%     of 1-bit (6dB) for most signals, as the exponent code is one larger
%     than it needs to be.
%
%     Conversely, what could be the reason for the choice? (1) Implementing
%     round-toward-zero requires adding the inverse of the sign to the
%     fractional component, and this logic might have be more expensive
%     than what Altera chose to implement. (2) Perhaps the FFT core
%     developer did not appreciate the impact of round-up when using block
%     floating-point format.
%
% 3. Burst I/O FFT Example Design with 16-bit data/twiddles
%
%    a) Input  data = 7FE0h (0111_1111_1110_0000b)
%       Output (data, exponent) = (03FFh, 5)
%
%    b) Input  data = 8000h (1000_0000_0000_0000b)
%       Output (data, exponent) = (FC00h, 5)
%
% This script shows the exponent code discrepancy between positive and
% negative impulse and constant inputs.
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
	pdfbase = sprintf('altera_burst_exp_codes_%d_%d', Bx, Bw);
else
	pdfbase = sprintf('radix4_bfp_exp_codes_%d_%d', Bx, Bw);
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
fft_analysis_exp_codes(fft_fixed, fft_double, Bx, N, pdfbase)
