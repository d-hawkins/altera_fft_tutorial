function fft_analysis_impulse(fft_fixed, fft_double, Bx, N, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_impulse.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Fixed-point Fast Fourier Transform (FFT) quantization analysis for impulse
% inputs.
%
% The script is for use with block floating-point fixed-point FFTs.
%
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------
%
% Figure numbers
fignum = 0;

% Extrema
x_max =  2^(Bx-1)-1;
%x_min = -2^(Bx-1);     % Asymmetric format
x_min = -x_max;       % Signed-symmetric format

% -----------------------------------------------------------------------------
% Impulse Analysis
% -----------------------------------------------------------------------------
%
fprintf('\n');
fprintf('Impulse Analysis\n');
fprintf('----------------\n');

fprintf(' * Create the impulses\n');

% Impulses
% - 16 phasors generated within a square (clockwise)
phi = pi/8;
r = round(x_max*tan(phi));
y = [
	 x_max   0
	 x_max  -r
	 x_max  x_min
	 r      x_min
	 0      x_min
	-r      x_min
	 x_min  x_min
	 x_min -r
	 x_min  0
	 x_min  r
	 x_min  x_max
	-r      x_max
	 0      x_max
	 r      x_max
	 x_max  x_max
	 x_max  r
];
M = length(y);

% Complex-valued
y = y(:,1) + 1j*y(:,2);

% Normalize to unity phasor
A = 2^(Bx-1);

% Polar plot
fignum = fignum+1;
figure(fignum)
clf
plot([1 1 -1 -1 1],[1 -1 -1 1 1],'k--')
hold on
for m = 1:M
	plot([0 real(y(m))]/A, [0 imag(y(m))]/A,'o-')
end
xlabel('Real')
ylabel('Imag')
axis(1.1*[-1 1 -1 1])
axis('square')
if strlength(pdfbase)
	pdfname = sprintf('%s_impulses.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

% FFT input vectors with impulse at sample x[1] (0-based indexing)
x = zeros(N,M);
for m = 1:M,
	% x[1] = impulse, for x[n], where n = 0, 1, ..., N-1
	x(2,m) = y(m);
	%
	% x[0] = impulse
%	x(1,m) = y(m);
	%
	% x[3] = impulse
%	x(4,m) = y(m);
	%
	% x[5] = impulse
%	x(6,m) = y(m);
end

fprintf(' * Calculate the FFTs\n');

% Spectra
X     = zeros(N,M);
X_exp = zeros(1,M);
Y     = zeros(N,M);
B_max = zeros(1,M);
for m = 1:M,
	X_est    = fft_fixed(x(:,m));
	X(:,m)   = X_est.data;
	X_exp(m) = X_est.exponent;
	Y(:,m)   = fft_double(x(:,m))/2^X_exp(m);

	% Maximum bit-widths (based on maximum magnitude)
	X_max = max([ max(abs(real(X(:,m)))) max(abs(imag(X(:,m)))) ]);
	B_max(m) = ceil(log2(X_max))+1;
end

% Magnitude bit-width required to represent the signed output
b = unique(B_max);
fprintf([' * Unique output bit-widths: [', ...
	repmat('%d ', 1, numel(b)-1), '%d]\n'],b)

% Print the unique exponent values
u = unique(X_exp);
fprintf([' * Unique exponent values: [', ...
	repmat('%d ', 1, numel(u)-1), '%d]\n'],u)

% Check the exponent width vs input data width
if ((Bx-u) < 6)
	fprintf('\nWARNING: an input bit-width of %d-bits is too small!\n\n', Bx)
end

fprintf(' * Plot the normalized impulse responses\n');

% Polar plot
fignum = fignum+1;
figure(fignum)
clf
%plot(X/2^(Bx-1))
plot(X(:,1)*2^X_exp(1)/A)
hold on
for m = 2:M
	plot(X(:,m)*2^X_exp(m)/A)
end
xlabel('Real')
ylabel('Imag')
axis(1.5*[-1 1 -1 1])
axis('square')
if strlength(pdfbase)
	pdfname = sprintf('%s_polar.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

% Real plot
n = [0:N-1];
fignum = fignum+1;
figure(fignum)
clf
subplot(2,1,1)
%plot(n,real(X/2^(Bx-1)))
plot(n, real(X(:,1))*2^X_exp(1)/A)
hold on
for m = 2:M
	plot(n, real(X(:,m))*2^X_exp(m)/A)
end
ylabel('Real')
xlabel('Frequency Channel')
axis([0 N-1 -1.5 1.5])
if strlength(pdfbase)
	pdfname = sprintf('%s_real.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

% Imag plot
fignum = fignum+1;
figure(fignum)
clf
subplot(2,1,1)
%plot(n,imag(X/2^(Bx-1)))
plot(n, imag(X(:,1))*2^X_exp(1)/A)
hold on
for m = 2:M
	plot(n, imag(X(:,m))*2^X_exp(m)/A)
end
ylabel('Imag')
xlabel('Frequency Channel')
axis([0 N-1 -1.5 1.5])
if strlength(pdfbase)
	pdfname = sprintf('%s_imag.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

fprintf(' * Plot the quantization error\n');

% Quantization error
Q = X-Y;

% Maximum error
Q_max = max([max(abs(real(Q))) max(abs(imag(Q)))]);
fprintf(' * Impulse quantization error %.2f\n', Q_max)

% Axis limits
Q_lim = 0.1*ceil((1.3*Q_max)/0.1);

% Quantization error polar plot
fignum = fignum+1;
figure(fignum)
clf
plot(Q,'.')
hold on
plot(Q_max*[1 1 -1 -1 1],Q_max*[1 -1 -1 1 1],'k--')
xlabel('Real')
ylabel('Imag')
axis(Q_lim*[-1 1 -1 1])
axis('square')
text(0.0,  1.1*Q_max, sprintf('%.2f', Q_max), ...
	'HorizontalAlignment','Center')
if strlength(pdfbase)
	pdfname = sprintf('%s_error_polar.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

% Individual error polar plots
for m = 1:M
	fignum = fignum+1;
	figure(fignum)
	clf
	plot((Q(:,m)),'.')
	hold on
	plot(Q_max*[1 1 -1 -1 1],Q_max*[1 -1 -1 1 1],'k--')
	xlabel('Real')
	ylabel('Imag')
	axis(Q_lim*[-1 1 -1 1])
	axis('square')
	text(0.0,  1.1*Q_max, sprintf('%.2f', Q_max), ...
		'HorizontalAlignment','Center')
	text(0.0, -1.1*Q_max, ...
		sprintf('Impulse (I, Q) = (%.2f, %.2f)', ...
		real(y(m))/A, imag(y(m))/A), ...
		'HorizontalAlignment','Center')
	if strlength(pdfbase)
		pdfname = sprintf('%s_error_polar_%d.pdf', pdfbase, m);
		exportgraphics(gca, pdfname)
	end
end

% Real plot
n = [0:N-1];
fignum = fignum+1;
figure(fignum)
clf
subplot(2,1,1)
plot(n,real(Q))
hold on
plot([0 N-1],+Q_max*[1 1],'k--')
plot([0 N-1],-Q_max*[1 1],'k--')
%text(50, 1.1*Q_max, sprintf('%.2f', Q_max))
text(50, 1.2*Q_max, sprintf('%.2f', Q_max))
ylabel('Real')
xlabel('Frequency Channel')
axis([0 N-1 -Q_lim Q_lim])
if strlength(pdfbase)
	pdfname = sprintf('%s_error_real.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

% Imag plot
fignum = fignum+1;
figure(fignum)
clf
subplot(2,1,1)
plot(n,imag(Q))
hold on
plot([0 N-1],+Q_max*[1 1],'k--')
plot([0 N-1],-Q_max*[1 1],'k--')
%text(50, 1.1*Q_max, sprintf('%.2f', Q_max))
text(50, 1.2*Q_max, sprintf('%.2f', Q_max))
ylabel('Imag')
xlabel('Frequency Channel')
axis([0 N-1 -Q_lim Q_lim])
if strlength(pdfbase)
	pdfname = sprintf('%s_error_imag.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

drawnow
