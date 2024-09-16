function fft_analysis_exp_codes(fft_fixed, fft_double, Bx, N, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_exp_codes.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Block Floating-Point Fast Fourier Transform (FFT) exponent code analysis
% using input impulses and constant values.
%
% -----------------------------------------------------------------------------

% Figure numbers
fignum = 0;

% -----------------------------------------------------------------------------
% Input impulses
% -----------------------------------------------------------------------------
%
% Number of impulses
M     = 10;

% Impulses and FFTs
x     = zeros(N,M);
y     = zeros(N,M);
X     = zeros(N,M);
X_exp = zeros(1,M);
Y     = zeros(N,M);
Y_exp = zeros(1,M);
for m = 1:M,
	% Impulse at x[0] (to generate constant spectra)
	%
	% Exponent = 6 at full-scale
%	x(1,m)   =  (2^(Bx-m)-1);
%	x(1,m)   =  1j*(2^(Bx-m)-1);
	x(1,m)   =  (2^(Bx-m)-1) + 1j*(2^(Bx-m)-1);
	%
	% Exponent = 5 at full-scale
%	y(1,m)   = -(2^(Bx-m));
%	y(1,m)   = -1j*(2^(Bx-m));
%	y(1,m)   = -(2^(Bx-m))   - 1j*(2^(Bx-m));
	%
	% Exponent = 5 at symmetric full-scale
%	y(1,m)   = -(2^(Bx-m)-1)   - 1j*(2^(Bx-m)-1);
	y(1,m)   = -x(1,m);

	% Fixed-point FFT
	X_est    = fft_fixed(x(:,m));
	X(:,m)   = X_est.data;
	X_exp(m) = X_est.exponent;

	% MATLAB FFT (error check)
	X_fft = fft(x(:,m));
	X_err = max(abs(X_est.data-X_fft/2^X_est.exponent));
	assert(X_err<1);

	% Fixed-point FFT
	Y_est    = fft_fixed(y(:,m));
	Y(:,m)   = Y_est.data;
	Y_exp(m) = Y_est.exponent;

	% MATLAB FFT (error check)
	Y_fft = fft(y(:,m));
	Y_err = max(abs(Y_est.data-Y_fft/2^Y_est.exponent));
	assert(Y_err<1);
end

fignum = fignum + 1;
figure(fignum)
clf
clear ph
ph(1) = plot(Bx-[0:M-1], X_exp, 'bx-');
hold on
ph(2) = plot(Bx-[0:M-1], Y_exp, 'ro-');
xlabel('Impulse Bitwidth')
ylabel('Exponent Code')
axis([Bx-M+0.5 Bx+0.5 -0.5 6.5])
legend(ph, ...
	'Positive Impulse',...
	'Negative Impulse',...
	'Location','NorthWest','AutoUpdate','Off')
if strlength(pdfbase)
	pdfname = sprintf('%s_impulse.pdf', pdfbase);
	pdfname = strrep(pdfname,'-','n');
	exportgraphics(gca, pdfname)
end

% -----------------------------------------------------------------------------
% Input constant
% -----------------------------------------------------------------------------
%
M     = 10;
x     = zeros(N,M);
y     = zeros(N,M);
X     = zeros(N,M);
X_exp = zeros(1,M);
Y     = zeros(N,M);
Y_exp = zeros(1,M);
for m = 1:M,
	% Constant values
	%
	% Exponent = 6 at full-scale
	x(:,m)   =  (2^(Bx-m)-1) + 1j*(2^(Bx-m)-1);
	%
	% Exponent = 5 at full-scale
%	y(:,m)   = -(2^(Bx-m))   - 1j*(2^(Bx-m));
	%
	% Exponent = 5 at symmetric full-scale
%	y(:,m)   = -(2^(Bx-m)-1)   - 1j*(2^(Bx-m)-1);
	y(:,m)   = -x(:,m);

	% Altera FFT
	X_est    = fft_fixed(x(:,m));
	X(:,m)   = X_est.data;
	X_exp(m) = X_est.exponent;

	X_fft = fft(x(:,m));
	X_err = max(abs(X_est.data-X_fft/2^X_est.exponent))/abs(X_est.data(1));
	assert(X_err<0.005); % less than 0.5-percent error

	% Altera FFT
	Y_est    = fft_fixed(y(:,m));
	Y(:,m)   = Y_est.data;
	Y_exp(m) = Y_est.exponent;

	Y_fft = fft(y(:,m));
	Y_err = max(abs(Y_est.data-Y_fft/2^Y_est.exponent))/abs(Y_est.data(1));
	assert(Y_err<0.005); % less than 0.5-percent error
end

fignum = fignum + 1;
figure(fignum)
clf
clear ph
ph(1) = plot(Bx-[0:M-1], X_exp, 'bx-');
hold on
ph(2) = plot(Bx-[0:M-1], Y_exp, 'ro-');
xlabel('Constant Bitwidth')
ylabel('Exponent Code')
axis([Bx-M+0.5 Bx+0.5 1.5 12.5])
legend(ph, ...
	'Positive Constant',...
	'Negative Constant',...
	'Location','NorthWest','AutoUpdate','Off')
if strlength(pdfbase)
	pdfname = sprintf('%s_constant.pdf', pdfbase);
	pdfname = strrep(pdfname,'-','n');
	exportgraphics(gca, pdfname)
end

