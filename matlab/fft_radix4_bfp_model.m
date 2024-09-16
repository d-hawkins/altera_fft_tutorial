% -----------------------------------------------------------------------------
% fft_radix4_bfp_model.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% FFT radix-4 model using double-precision calculations.
%
% -----------------------------------------------------------------------------
function y = fft_radix4_bfp_model(x, Bx, Bw, Bg)

	% Require Nx1 input (and return Nx1)
	assert(size(x,1)~=1)

	% Default guard bits
	if (nargin < 4)
		Bg = 3;
	end

	% Rounding
	% * "tozero", "fromzero", "odd", "even"
	% * "tozero" and "fromzero" create quantization noise bias
	% * "odd" and "even" create uniform quantization noise
	% * "minusinf" or "plusinf" for scale_round causes a DC spike
	% * "truncate" truncate (bad!)
	% * "add_half" add half and truncate (not as bad)
	scale_round   = "odd";
	product_round = "odd";

	% Length
	N = length(x);

	% Number of stages
	S = log(N)/log(4);

	% Twiddles table
	n = [0:N-1]';
%	W = exp(-j*2*pi/N*n);
	W = round((2^(Bw-1)-1)*exp(-j*2*pi/N*n))/2^(Bw-1);
%	W = round(2^(Bw-1)*exp(-j*2*pi/N*n))/2^(Bw-1);

	% In-place initialize
	X = x;

	% Process all stages
	X_exp = 0;
	for s = 0:S-1,

		% Number of groups
		% * increases by 4 for each stage
		G = 4^s;

		% Number of butterflies per group
		% * decreases by 4 for each stage
		B = 4^(S-1-s);

		% Twiddles for a single group
		% - Scaled by G to the N-point roots-of-unity phasor
		n = [0:B-1];
		k = [0:3]';
		kn = mod(G*k*n,N);

%		fprintf('S = %d, G = %d, B = %d\n', s, G, B)

		% Block floating-point scaling
		% * Scale before processing each stage
		% * This is equivalent to a hardware implementation that
		%   determines the required scale value on data load and
		%   then applies it before the first pass, and then during
		%   subsequent passes calculates the next required scaling
		scale = bfp_scale(Bx, Bg, X);
		if (scale)
%			fprintf(' * scale = %d\n', scale)

			% Divide-by-2, 4, or 8
			switch scale_round
				case "truncate"
					X = floor(X/2^scale);

				case "add_half"
					% Add 0.5 and truncate - creates a DC spike
					X = floor(X/2^scale + 0.5*(1+1j)*ones(N,1));

				otherwise
					X = round(X/2^scale, TieBreaker=scale_round);
			end

			% Increment the exponent
			X_exp = X_exp + scale;
		end

		% Loop over the groups (columns) of radix-4 butterflies
		for g = 0:G-1,
			% Bufferflies
			for b = 0:B-1,
%				fprintf('  * g = %2d, b = %2d\n', g, b)

				% In-place data indices for a butterfly group
				% * 4 indices within the columns of a 2D breakdown
				m = 4*B*g + [0:3]'*B + b;

%				fprintf('  * indices = (%d, %d, %d, %d)\n', ...
%					m(1), m(2), m(3), m(4))

				X_4 = X(m+1);
				X_4 = bfy_radix4(X_4);

				% Twiddles
				if (s < (S-1))
%					fprintf('  * twiddles = (%d, %d, %d, %d)\n', ...
%						kn(1,b+1), kn(2,b+1), kn(3,b+1), kn(4,b+1))
					X_4 = X_4.*W(kn(:,b+1)+1);
				end

				switch product_round
					case "truncate"
						% Causes a DC spike and passband spikes
						X_4 = floor(X_4);

					case "add_half"
						% Add 0.5 and truncate
						X_4 = floor(X_4 + 0.5*(1+1j)*ones(4,1));

					otherwise
						X_4 = round(X_4, TieBreaker=product_round);
				end

				% In-place write-back
				X(m+1) = X_4;

			end
		end
	end

	% Digit reverse
	X = digitrevorder(X,4);

	% Return structure
	y.data     = X;
	y.exponent = X_exp;
end

% -----------------------------------------------------------------------------
% Radix-4 (4-point DFT) butterfly
% -----------------------------------------------------------------------------
%
function y = bfy_radix4(x)
	y = [1 1 1 1; 1 -j -1 j; 1 -1 1 -1; 1 j -1 -j]*x;
end

% -----------------------------------------------------------------------------
% Block Floating-Point (BFP) scaling
% -----------------------------------------------------------------------------
%
function scale = bfp_scale(Bx, Bg, X)
	% Default (no scaling)
	scale = 0;

	% Maximum integer value
	X_max = max([ max(abs(real(X))) max(abs(real(X))) ]);

	% There should be Bg = 3 spare magnitude bits
	spare_bits = Bx - 1 - Bg - ceil(log2(X_max));
	if (spare_bits < 0)
		scale = -spare_bits;
	end
end
