function fft_analysis_npr(fft_fixed, fft_double, Bx, N, M, LF_dB, fft_fixed_name, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_npr.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Fixed-point Fast Fourier Transform (FFT) quantization analysis for
% complex-valued Gaussian notch-filtered noise for a range of input
% loading factors.
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
% Number of coefficients in the filter
Nh = 512;

% Filter coefficients response window
wh = kaiser(Nh,25).';

% Passband bandwidth in channels
Bh = 2*round(0.9*Nh/2);

% Stopband bandwidth in channels
Sh = 2*round(0.1*Nh/2);

% Slope over the band
SdB = [0  0];
%SdB = [0 20];

% Baseband (ideal) frequency response (in dB)
% * the zero channel is at index (Nh/2+1)
H_dB = -200*ones(1,Nh); % -200dB
%
% Passband
H_dB(Nh/2+1+[-Bh/2:Bh/2]) = SdB(1) + [0:Bh]/Bh*(SdB(2)-SdB(1));
%
% Stopband
H_dB(Nh/2+1+[-Sh/2:Sh/2]) = -200;

% Linear response
H = 10.^(H_dB/20);

% Windowed coefficients response (complex valued)
h = wh.*fftshift(ifft(fftshift(H)));

% Power spectra
H = fftshift(fft(fftshift(h)));
h_ig = sqrt(mean(abs(H).^2));
H_dB = 10*log10(abs(H).^2 + 10^(-20));

% -----------------------------------------------------------------------------
% Gaussian Notch-Filtered Noise NPR Analysis
% -----------------------------------------------------------------------------
%
fprintf('\n');
fprintf('Gaussian Notch-Filtered Noise NPR Analysis\n');
fprintf('------------------------------------------\n');

fprintf(' * Create the complex-valued noise samples\n');

% Random number generator
rng(1234,"twister");

% Total number of complex-valued samples to filter
Ns = N*M + Nh;

% Complex-valued noise
x_double = randn(1, Ns) + 1j*randn(1, Ns);

% Filter the noise
x_double = filter(h,1,x_double);

% Discard the transient
x_double = x_double(Nh + [1:N*M]);

% RMS
x_std = std(x_double);

% Reshape (for windowing)
x_double = reshape(x_double,N,M);

% Reference double-precision data used each loop
x_reference = x_double;

% -----------------------------------------------------------------------------
% Loop over loading factors
% -----------------------------------------------------------------------------
%
% Loading factors list
LF = 10.^(LF_dB./20);

% Loop over loading factors
K = length(LF_dB);
X_npr_dB  = zeros(1,K);
X_exp_avg = zeros(3,K);
Y_npr_dB  = zeros(1,K);
for k = 1:K,

	fprintf(' * Loading factor = %.2f dB\n', LF_dB(k));

	% Normalize and scale to the desired loading factor
	x_double = LF(k)*2^(Bx-1)*x_reference/x_std;

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
%		Y(:,m) = fft_double(x_double(:,m));

		% Double-precision FFT (of the fixed-point input)
		% * Use this to calculate the fixed-point difference noise
		% * Required for NPR calculation at the band center
		%   (see Y_pass_dB, Y_stop_dB and Y_npr_dB)
		Y(:,m) = fft_double(x_fixed(:,m));

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

	X_exp_avg(1,k) = max(X_exp);
	X_exp_avg(2,k) = mean(X_exp);
	X_exp_avg(3,k) = min(X_exp);

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

	% Noise power ratio
	N_avg = 2*round(Sh/Nh*N/4);
	N_pass = (N/2+1)+N/4+[-N_avg/2:N_avg/2];
	N_stop = (N/2+1)+[-N_avg/2:N_avg/2];

	% Fixed point FFT
	X_pass_dB = 10*log10(mean(Rxx(N_pass)));
	X_stop_dB = 10*log10(mean(Rxx(N_stop)));
	X_npr_dB(k) = X_pass_dB-X_stop_dB;

	% Floating point FFT
	Y_pass_dB = 10*log10(mean(Ryy(N_pass)));
	Y_stop_dB = 10*log10(mean(Ryy(N_stop)));
	Y_npr_dB(k) = Y_pass_dB-Y_stop_dB;

	fprintf(' * Noise Power Ratio (NPR)\n')
	fprintf('    - Block floating-point FFT NPR = %.2f dB\n', X_npr_dB(k))
	fprintf('    - Double-precision     FFT NPR = %.2f dB\n', Y_npr_dB(k))

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
	ph(1) = plot(f, Ryy_dB, 'g');
	ph(2) = plot(f, Rxx_dB, 'b');
	ph(3) = plot(f, Rqq_dB, 'r');
	plot(f([1 end]), Rxx_mean_dB*[1 1], 'b--')
	plot(f([1 end]), Q_dB*[1 1], 'm--')
	xlabel('Frequency Channel')
	ylabel('Loading Factor (dB)')
	y_min = 10*floor(Q_dB/10);
	%y_min = 10*floor(min(Rqq_dB)/10);
	%y_min = -80;
	axis([f([1 end]) y_min 10])
	legend(ph,'MATLAB FFT', fft_fixed_name, 'Difference', ...
		'Location','NorthEast','AutoUpdate','Off')

	% Input data loading factor
	text(0, Rxx_mean_dB+3, sprintf('%.2f dB', Rxx_mean_dB),'Color','b', ...
		'HorizontalAlignment','Center')

	% Plot the maximum value of the quantization noise
	text(-N/4, Rqq_max_dB+3, sprintf('%.2f dB', Rqq_max_dB),'Color','r', ...
		'HorizontalAlignment','Center')

	% Input quantization noise floor
	text(0, Q_dB+3, sprintf('%.2f dB', Q_dB),'Color','m', ...
		'HorizontalAlignment','Center')

	%axis([f([1 end]) -112 10])

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

	% Re-use the same figure numbers for each loop
	if (k < K)
		fignum = fignum-4;
	end
end

% Noise Power Ratio
fignum = fignum+1;
figure(fignum)
clf
clear ph
% MATLAB FFT NPR
ph(1) = plot(LF_dB, Y_npr_dB, 'r.-');
hold on
%
% Fixed-point FFT NPR
% * The Radix-4 BFP NPR is slightly lower for LF < -40dB due to
%   the rounding of the twiddle factors in the model
% * The Altera model is lower by about 2-bits
ph(2) = plot(LF_dB, X_npr_dB, 'b.-');
%
% Maximum noise power ratio for the MATLAB FFT
Y_npr_max_dB = max(Y_npr_dB);
plot(LF_dB([1 end]),Y_npr_max_dB*[1 1],'r--')
text(-35, Y_npr_max_dB-2.5, sprintf('%.2f dB', Y_npr_max_dB),'Color','r')
%
% Maximum noise power ratio for the fixed-point FFT
X_npr_max_dB = max(X_npr_dB);
plot(LF_dB([1 end]),X_npr_max_dB*[1 1],'b--')
text(-15, X_npr_max_dB+3, sprintf('%.2f dB', X_npr_max_dB),'Color','b')
%
xlabel('Loading Factor (dB)')
ylabel('Noise Power Ratio (dB)')
y_max = 10*ceil(max(Y_npr_dB)/10+0.5);
axis([LF_dB([1 end]) 0 y_max])
legend(ph,'MATLAB FFT', fft_fixed_name, ...
	'Location','SouthWest','AutoUpdate','Off')
if strlength(pdfbase)
	pdfname = sprintf('%s.pdf', pdfbase);
	pdfname = strrep(pdfname,'-','n');
	exportgraphics(gca, pdfname)
end

% Exponent code
% * The slope should change at about 6dB/bit
fignum = fignum+1;
figure(fignum)
clf
clear ph
%
% Plot a 6.02dB/bit slope
m = find(X_exp_avg(2,:) > 3);
m = m([1:end-1]); % Discard the last sample
offset = mean(X_exp_avg(2,m)-LF_dB(m)/6.02);
plot(LF_dB, LF_dB/6.02+offset, 'k--');
hold on
ph(1) = plot(LF_dB, X_exp_avg(1,:), 'rv-');
ph(3) = plot(LF_dB, X_exp_avg(3,:), 'g^-');
% Plot the mean last, so that it overlays min and max
ph(2) = plot(LF_dB, X_exp_avg(2,:), 'bo-','LineWidth',1.0,'MarkerFaceColor','b');
xlabel('Loading Factor (dB)')
ylabel('Exponent Code')
legend(ph,'Max','Mean','Min', ...
	'Location','NorthWest','AutoUpdate','Off')
y_max = ceil(max(X_exp_avg(2,:)));
%axis([ LF_dB([1 end])+[-1 1] -0.1 y_max])
axis([LF_dB([1 end]) 0 10])
if strlength(pdfbase)
	pdfname = sprintf('%s_exp_codes.pdf', pdfbase);
	pdfname = strrep(pdfname,'-','n');
	exportgraphics(gca, pdfname)
end

%keyboard

