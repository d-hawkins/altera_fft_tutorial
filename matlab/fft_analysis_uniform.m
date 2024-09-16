function fft_analysis_uniform(fft_fixed, fft_double, Bx, N, M, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_uniform.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Fixed-point Fast Fourier Transform (FFT) quantization analysis for
% complex-valued uniform noise inputs.
%
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------
%
% Figure numbers
fignum = 0;

% -----------------------------------------------------------------------------
% Uniform Noise Analysis
% -----------------------------------------------------------------------------
%
fprintf('\n');
fprintf('Uniform Noise Analysis\n');
fprintf('----------------------\n');

fprintf(' * Create the complex-valued noise samples\n');

% Reset the random number generator
rng(1234,'twister');

% Real-valued random data in double-precision format
x_double = 2^(Bx-1)*(rand(2*M*N,1)-0.5);

% Quantized
x_fixed = round(x_double);

% Saturated
x_max = 2^(Bx-1)-1;
m = find(x_fixed > x_max);
if ~isempty(m)
	fprintf('Saturate %d +ve samples\n', length(m))
	x_fixed(m) = x_max;
end
m = find(x_fixed < -x_max);
if ~isempty(m)
	fprintf('Saturate %d -ve samples\n', length(m))
	x_fixed(m) = -x_max;
end

% Complex-valued double-precision
x_double = reshape(x_double(1:2:2*M*N),N,M) + ...
	1j*reshape(x_double(2:2:2*M*N),N,M);

% Complex-valued fixed-point
x_fixed = reshape(x_fixed(1:2:2*M*N),N,M) + ...
	1j*reshape(x_fixed(2:2:2*M*N),N,M);

% Quantization error polar plot
fignum = fignum+1;
figure(fignum)
clf
plot((x_fixed-x_double),'.')
xlabel('Real')
ylabel('Imag')
axis(0.6*[-1 1 -1 1])
axis('square')
title('Uniform Noise Quantization Error')

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
fprintf(' * Uniform noise quantization noise floor %.2f dB\n', Rqq_max_dB)

% Frequency channels in frequency-domain order
f = [-N/2:N/2-1];

% Quantization noise floor for complex-valued data
Q_dB = -(6.04*Bx+1.76);

% Spectra
fignum = fignum+1;
figure(fignum)
clf
clear ph
plot(f([1 end]), [0 0], 'k--')
hold on
ph(1) = plot(f, Rxx_dB, 'b');
ph(2) = plot(f, Rqq_dB, 'r');
plot(f([1 end]), Q_dB*[1 1], 'm--')
xlabel('Frequency Channel')
ylabel('Loading Factor (dB)')
y_min = 10*floor(Q_dB/10);
%y_min = 10*floor(min(Rqq_dB)/10);
%y_min = -80;
axis([f([1 end]) y_min 10])
legend(ph,'Input Uniform Noise', 'Quantization Noise', ...
	'Location','East','AutoUpdate','Off')

% Input data loading factor
text(0, Rxx_mean_dB+3, sprintf('%.2f dB', Rxx_mean_dB),'Color','b', ...
	'HorizontalAlignment','Center')

% Plot the maximum value of the quantization noise
text(0, Rqq_max_dB+3, sprintf('%.2f dB', Rqq_max_dB),'Color','r', ...
	'HorizontalAlignment','Center')

% Input quantization noise floor
text(0, Q_dB+3, sprintf('%.2f dB', Q_dB),'Color','m', ...
	'HorizontalAlignment','Center')

%axis([f([1 end]) -112 10])
if strlength(pdfbase)
	pdfname = sprintf('%s_spectra.pdf', pdfbase);
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

