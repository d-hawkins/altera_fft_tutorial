function fft_analysis_exponential(fft_fixed, fft_double, Bx, N, pdfbase)
% -----------------------------------------------------------------------------
% fft_analysis_exponential.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% Fixed-point Fast Fourier Transform (FFT) quantization analysis for
% complex-valued exponent inputs.
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

% Exponent amplitude (full-scale)
A =  2^(Bx-1)-1;

% -----------------------------------------------------------------------------
% Impulse Analysis
% -----------------------------------------------------------------------------
%
fprintf('\n');
fprintf('Complex-valued Exponent Analysis\n');
fprintf('--------------------------------\n');

fprintf(' * Create the random frequency and random phase exponent signals\n');

% Number of complex-valued exponential inputs
M = 128;

% Reset the random number generator
rng(12345,'twister');

% Random center frequencies (uniform distribution)
f = round(rand(1,M)*N);

% Random phase (uniform distribution)
ph = 2*pi*rand(1,M);

% Or: Uniform phase steps
%ph = 2*pi/M*[0:M-1];

% Sample index
n = [0:N-1];

% FFT inputs
M = length(f);
x = zeros(N,M);
for m = 1:M,
	x(:,m) = round(A*exp(j*(2*pi*n/N*f(m)+ph(m))));
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

fprintf(' * Plot the normalized impulse responses\n');

% Polar plot
fignum = fignum+1;
figure(fignum)
clf
%plot(X/2^(Bx-1))
X_norm = X(:,1)*2^X_exp(1)/(A*N);
plot(real(X_norm),imag(X_norm))
hold on
for m = 2:M
	X_norm = X(:,m)*2^X_exp(m)/(A*N);
	plot(real(X_norm),imag(X_norm))
end
xlabel('Real')
ylabel('Imag')
axis(1.1*[-1 1 -1 1])
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
%plot(n,real(X/2^(Bx-1)))
plot(n, real(X(:,1))*2^X_exp(1)/(A*N))
hold on
for m = 2:M
	plot(n, real(X(:,m))*2^X_exp(m)/(A*N))
end
ylabel('Real')
xlabel('Frequency Channel')
axis([0 N-1 -1.1 1.1])

% Imag plot
fignum = fignum+1;
figure(fignum)
clf
%plot(n,imag(X/2^(Bx-1)))
plot(n, imag(X(:,1))*2^X_exp(1)/(A*N))
hold on
for m = 2:M
	plot(n, imag(X(:,m))*2^X_exp(m)/(A*N))
end
ylabel('Imag')
xlabel('Frequency Channel')
axis([0 N-1 -1.1 1.1])

% Mag plot
fignum = fignum+1;
figure(fignum)
clf
plot(n, abs(X(:,1))*2^X_exp(1)/(A*N))
hold on
for m = 2:M
	plot(n, abs(X(:,m))*2^X_exp(m)/(A*N))
end
ylabel('Magnitude')
xlabel('Frequency Channel')
axis([0 N-1 -0.1 1.1])
if strlength(pdfbase)
	pdfname = sprintf('%s_channels.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

fprintf(' * Plot the quantization error\n');

% Quantization error
Q = X-Y;

% Maximum error
Q_max = max([max(abs(real(Q))) max(abs(imag(Q)))]);
fprintf(' * Quantization error %.2f\n', Q_max)

% Axis limits
Q_lim = 0.5*ceil((1.2*Q_max)/0.5);

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
text(-0.1, 1.1*Q_max, sprintf('%.2f', Q_max))
if strlength(pdfbase)
	pdfname = sprintf('%s_error.pdf', pdfbase);
	exportgraphics(gca, pdfname)
end

if (0)
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
		text(-0.1, 1.1*Q_max, sprintf('%.2f', Q_max))
	end
end

% Real plot
n = [0:N-1];
fignum = fignum+1;
figure(fignum)
clf
plot(n,real(Q))
hold on
plot([0 N-1],+Q_max*[1 1],'k--')
plot([0 N-1],-Q_max*[1 1],'k--')
text(50, 1.1*Q_max, sprintf('%.2f', Q_max))
ylabel('Real')
xlabel('Frequency Channel')
axis([0 N-1 -Q_lim Q_lim])

% Imag plot
fignum = fignum+1;
figure(fignum)
clf
plot(n,imag(Q))
hold on
plot([0 N-1],+Q_max*[1 1],'k--')
plot([0 N-1],-Q_max*[1 1],'k--')
text(50, 1.1*Q_max, sprintf('%.2f', Q_max))
ylabel('Imag')
xlabel('Frequency Channel')
axis([0 N-1 -Q_lim Q_lim])

drawnow
