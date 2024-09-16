function fft_analysis_bandpass(fft_fixed, fft_double, Bx, N, M, LF_dB, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_bandpass.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Fixed-point Fast Fourier Transform (FFT) quantization analysis for
% complex-valued Gaussian bandpass noise.
%
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------
%
% Figure numbers
fignum = 0;

% -----------------------------------------------------------------------------
% FFT window function
% -----------------------------------------------------------------------------
%
% The frequency channel response of the FFT is sinc-like. This sinc-like
% response has such high sidelobes, that the stopband of the filtered input
% noise signal gets filled in, so you cannot see the quantization noise floor.
% This m-file windows the floating-point data and then quantizes it. This
% reduces the sidelobe response, but leaves the quantization noise floor where
% it is expected.
%

% Window function
w = kaiser(N+1,25);
w = w(1:N);      % DFT cyclic symmetry

% Power spectrum normalization
w_cg    = mean(w);          % Coherent gain
w_ig    = sqrt(mean(w.^2)); % Incoherent gain
w_cg_dB = 20*log10(w_cg);
w_ig_dB = 20*log10(w_ig);

% -----------------------------------------------------------------------------
% Complex-valued baseband filter
% -----------------------------------------------------------------------------
%
% The noise signal is filtered with a complex-valued baseband signal defined
% by a filter with a different negative and positive frequency response.
% * Positive: bandpass with a linear ramp over the passband
% * Negative: bandpass with a semi-circular passband
%

% Number of coefficients in the filter
Nh = 512;

% Filter coefficients response window
wh = kaiser(Nh,25).';

% Bandwidth in channels
%Bh = Nh/4;
Bh = 2*round(Nh/8*1.2);

% Filter#1: passband ramp
% -----------------------
%
% Slope over the band
%SdB = [0  0];
SdB = [0 20];
%
% Baseband (ideal) frequency response (in dB)
% * the zero channel is at index (Nh/2+1)
H1_dB = -200*ones(1,Nh); % -200dB
H1_dB(Nh/2+1+[-Bh/2:Bh/2]) = SdB(1) + [0:Bh]/Bh*(SdB(2)-SdB(1));

% Filter#2: passband semi-circle
% ------------------------------
%
% Adjust the offset to align the semi-circle passband
% with the sloped response
offset = -58;
%offset = -48; % good for Bh = Nh/4
H2_dB = -200*ones(1,Nh); % -200dB
H2_dB(Nh/2+1+[-Bh/2:Bh/2]) = sqrt( (Bh/2)^2 - [-Bh/2:Bh/2].^2) + offset;

% Filter complex-valued baseband responses
% ----------------------------------------
%
% Linear responses
H1 = 10.^(H1_dB/20);
H2 = 10.^(H2_dB/20);

% Windowed coefficients response (complex valued)
h1 = wh.*fftshift(ifft(fftshift(H1)));
h2 = wh.*fftshift(ifft(fftshift(H2)));

% Power spectra
H1 = fftshift(fft(fftshift(h1)));
H2 = fftshift(fft(fftshift(h2)));
h1_ig = sqrt(mean(abs(H1).^2));
h2_ig = sqrt(mean(abs(H2).^2));

H1_dB = 10*log10(abs(H1).^2 + 10^(-20));
H2_dB = 10*log10(abs(H2).^2 + 10^(-20));

% Filter total response
% ---------------------
%
% Modulate each filter to +/- Nh/4
nh = [-Nh/2:Nh/2-1];
h  = h1.*exp(j*2*pi*nh/4) + h2.*exp(-j*2*pi*nh/4);

H = fftshift(fft(fftshift(h)));
h_ig = sqrt(mean(abs(H).^2));
H_dB = 10*log10(abs(H).^2 + 10^(-20));

% -----------------------------------------------------------------------------
% Gaussian Bandpass Noise Analysis
% -----------------------------------------------------------------------------
%
fprintf('\n');
fprintf('Gaussian Bandpass Noise Analysis\n');
fprintf('--------------------------------\n');

fprintf(' * Create the complex-valued noise samples\n');

% Random number generator
rng(1234,"twister");

% Loading factor
LF = 10^(LF_dB./20);

% Total number of complex-valued samples to filter
Ns = N*M + Nh;

% Complex-valued noise
x_double = randn(1, Ns) + 1j*randn(1, Ns);

% Filter the noise
x_double = filter(h,1,x_double);

% Discard the transient
x_double = x_double(Nh + [1:N*M]);

% Normalize and scale to the desired loading factor
x_std    = std(x_double);
x_double = LF*2^(Bx-1)*x_double/x_std;

% Convert to N x M
x_double = reshape(x_double,N,M);

% Apply the FFT window function to each of the M estimates
for m = 1:M,
	x_double(:,m) = w.*x_double(:,m)/w_ig;
end

% Quantize
x_fixed = round(x_double);

% Saturate
x_max = 2^(Bx-1)-1;
x_fixed_re = real(x_fixed);
x_fixed_im = imag(x_fixed);
for m = 1:M,
	p = find(x_fixed_re(:,m) > x_max);
	if ~isempty(p)
		fprintf('Saturate %d +ve real-samples\n', length(p))
		x_fixed_re(p,m) = x_max;
	end
	p = find(x_fixed_re(:,m) < -x_max);
	if ~isempty(p)
		fprintf('Saturate %d -ve real-samples\n', length(p))
		x_fixed_re(p,m) = -x_max;
	end
	p = find(x_fixed_im(:,m) > x_max);
	if ~isempty(p)
		fprintf('Saturate %d +ve imag-samples\n', length(p))
		x_fixed_im(p,m) = x_max;
	end
	p = find(x_fixed_im(:,m) < -x_max);
	if ~isempty(p)
		fprintf('Saturate %d -ve imag-samples\n', length(p))
		x_fixed_im(p,m) = -x_max;
	end
end
x_fixed = x_fixed_re + 1j*x_fixed_im;
clear x_fixed_re x_fixed_im

% Quantization error polar plot
fignum = fignum+1;
figure(fignum)
clf
plot((x_fixed-x_double),'.')
xlabel('Real')
ylabel('Imag')
axis(0.6*[-1 1 -1 1])
axis('square')
title('Gaussian Noise Quantization Error')

fprintf(' * Calculate the FFTs\n');


% Spectra
X     = zeros(N,M);
X_exp = zeros(1,M);
Y     = zeros(N,M);
B_max = zeros(1,M);
for m = 1:M,
	% Fixed-point FFT
	X_est    = fft_fixed(x_fixed(:,m));
	X_exp(m) = X_est.exponent;
	X(:,m)   = X_est.data*2^X_exp(m);

	% Double-precision FFT (of the floating-point input)
	% * Use this to calculate the total quantization noise
	% * Aligns with the noise floor of the BFP spectra
	Y(:,m) = fft_double(x_double(:,m));

	% Double-precision FFT (of the fixed-point input)
	% * Use this to calculate the fixed-point difference noise
%	Y(:,m) = fft_double(x_fixed(:,m));

	% Maximum bit-widths (based on maximum magnitude)
	X_max = max([ max(abs(real(X_est.data))) max(abs(imag(X_est.data))) ]);
	B_max(m) = ceil(log2(X_max))+1;

	% Progress indicator
	if (mod(m,20)==0)
		fprintf('.')
	end
	if (mod(m,500)==0)
		fprintf(' %d\n', m)
	end
end
fprintf('\n')

% Magnitude bit-width required to represent the signed output
b = unique(B_max);
fprintf([' * Unique output bit-widths: [', ...
	repmat('%d ', 1, numel(b)-1), '%d]\n'],b)

% Print the unique exponent values
u = unique(X_exp);
fprintf([' * Unique exponent values: [', ...
	repmat('%d ', 1, numel(u)-1), '%d]\n'],u)

% Quantization error
Q = X-Y;

fprintf(' * Calculate the power spectral estimates\n');

% Spectra
Rxx = zeros(N,1);
Ryy = zeros(N,1);
Rqq = zeros(N,1);
for m = 1:M,
	Rxx = Rxx + abs(X(:,m)).^2;
	Ryy = Ryy + abs(Y(:,m)).^2;
	Rqq = Rqq + abs(Q(:,m)).^2;
end
% Normalize
% * M for the average
% * (2^(Bx-1))^2 to convert to fractional integer
% * N for the FFT power gain
Rxx = Rxx/(M*2^(2*(Bx-1))*N);
Ryy = Ryy/(M*2^(2*(Bx-1))*N);
Rqq = Rqq/(M*2^(2*(Bx-1))*N);

% Convert spectra to frequency-domain order
Rxx = fftshift(Rxx);
Ryy = fftshift(Ryy);
Rqq = fftshift(Rqq);

% dB format
Rxx_dB = 10*log10(Rxx+10^(-100));
Ryy_dB = 10*log10(Ryy+10^(-100));
Rqq_dB = 10*log10(Rqq+10^(-100));

% Average noise power
Rxx_mean_dB = 10*log10(mean(Rxx));
Rqq_mean_dB = 10*log10(mean(Rqq));
Rqq_max_dB  = 10*log10(max(Rqq));
fprintf(' * Gaussian noise quantization noise floor %.2f dB\n', Rqq_max_dB)

% The index is needed for annotating the text in the right location
[~,Rqq_max_n] = max(Rqq);

% Frequency channels in frequency-domain order
f = [-N/2:N/2-1];

% Quantization noise floor for complex-valued data
Q_dB = -(6.04*Bx+1.76);

% Quantization noise floor for Radix-4 BFP bit-width
%  * The BFP signal-to-quantization noise is 3-bits less
%  * The noise floor is about 3dB higher due to the cascade
%    of scaling and product quantization. You can see this 3dB
%    offset disappear when the input signal loading factor is
%    small enough that the exponent code is 0, eg.,-60dB.
%  * BFP SNR is nominally (Bx-3.5)-bits
%B_dB = -(6.04*(Bx-3.5)+1.76);

% Spectra
fignum = fignum+1;
figure(fignum)
clf
clear ph
plot(f([1 end]), [0 0], 'k--')
hold on
ph(1) = plot(f, Rxx_dB, 'b');
ph(2) = plot(f, Rqq_dB, 'r');
plot(f([1 end]), Rxx_mean_dB*[1 1], 'b--')
plot(f([1 end]), Q_dB*[1 1], 'm--')
%plot(f([1 end]), B_dB*[1 1], 'c--')
xlabel('Frequency Channel')
ylabel('Loading Factor (dB)')
y_min = 10*floor(Q_dB/10);
%y_min = 10*floor(min(Rqq_dB)/10);
%y_min = -80;
axis([f([1 end]) y_min 10])
legend(ph,'Bandpass Gaussian Noise', 'Quantization Noise', ...
	'Location','North','AutoUpdate','Off')

% Input data loading factor
text(0, Rxx_mean_dB+3, sprintf('%.2f dB', Rxx_mean_dB),'Color','b', ...
	'HorizontalAlignment','Center')

if (abs(Rqq_max_dB-Q_dB)>10)
	% Plot the maximum value of the quantization noise
	text(f(Rqq_max_n), Rqq_max_dB+3, sprintf('%.2f dB', Rqq_max_dB),'Color','r', ...
		'HorizontalAlignment','Center')

	% Input quantization noise floor
	text(0, Q_dB+3, sprintf('%.2f dB', Q_dB),'Color','m', ...
		'HorizontalAlignment','Center')
else
	y_dB = max([Rqq_max_dB Q_dB]);

	% Plot the maximum value of the quantization noise
	text(-N/4, y_dB+3, sprintf('%.2f dB', Rqq_max_dB),'Color','r', ...
		'HorizontalAlignment','Center')

	% Input quantization noise floor
	text(+N/4, y_dB+3, sprintf('%.2f dB', Q_dB),'Color','m', ...
		'HorizontalAlignment','Center')
end

%axis([f([1 end]) -112 10])

% Noise power ratio
npr_dB = Rxx_mean_dB - Rqq_max_dB;
text(0, Rxx_mean_dB-5, sprintf('NPR = %.2f dB', npr_dB),'Color','k', ...
	'HorizontalAlignment','Center')

%axis([f([1 end]) -112 10])
if strlength(pdfbase)
	pdfname = sprintf('%s_%ddB_spectra.pdf', pdfbase, LF_dB);
	pdfname = strrep(pdfname,'-','n');
	exportgraphics(gca, pdfname)
end

% Error histogram
% * Uniform quantization noise RMS is 1/sqrt(12)
% * Real and imag spectra RMS should be sqrt(N) larger
%
% Normalized components
Q_re = real(Q(:))/sqrt(N/12);
Q_im = imag(Q(:))/sqrt(N/12);
Q_re_std = std(Q_re);
Q_im_std = std(Q_im);
fprintf(' * Normalized quantization noise RMS (re, im) = (%.3f, %.3f)\n', ...
	Q_re_std, Q_im_std);

% Histogram bin centers
c_max = 10*ceil(3.5*max([Q_re_std Q_im_std])/10);
c = linspace(-c_max,c_max,101);
Q_re_h = hist(Q_re,c);
Q_im_h = hist(Q_im,c);

% Maximum (for zooming both re and im histograms)
h_max = max([max(Q_re_h) max(Q_im_h)]);
h_max_log10 = floor(log10(h_max));
h_max = ceil(h_max/10^h_max_log10)*10^h_max_log10;

fignum = fignum+1;
figure(fignum)
clf
bar(c,Q_re_h,1,'FaceColor','b')
xlabel('Normalized Error Real-Part')
axis([-c_max c_max 0 h_max])

fignum = fignum+1;
figure(fignum)
clf
bar(c,Q_im_h,1,'FaceColor','b')
xlabel('Normalized Error Imag-Part')
axis([-c_max c_max 0 h_max])

drawnow

%keyboard

