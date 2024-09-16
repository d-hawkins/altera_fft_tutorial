% -----------------------------------------------------------------------------
% fft_radix4_model.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% FFT radix-4 model using double-precision calculations.
%
% -----------------------------------------------------------------------------
function X = fft_radix4_model(x,Bw)

	% Require Nx1 input (and return Nx1)
	assert(size(x,1)~=1)

	% Length
	N = length(x);

	% Number of stages
	S = log(N)/log(4);

	% Twiddles table
	n = [0:N-1]';
%	W = exp(-j*2*pi/N*n);
	W = (2^(Bw-1)-1)*exp(-j*2*pi/N*n)/2^(Bw-1);

	% In-place initialize
	X = x;

	% Process all stages
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

		fprintf('S = %d, G = %d, B = %d\n', s, G, B)

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

				% In-place write-back
				X(m+1) = X_4;

			end
		end
	end

	% Digit reverse
	X = digitrevorder(X,4);
end

% -----------------------------------------------------------------------------
% Radix-4 (4-point DFT) butterfly
% -----------------------------------------------------------------------------
%
function y = bfy_radix4(x)
	y = [1 1 1 1; 1 -j -1 j; 1 -1 1 -1; 1 j -1 -j]*x;
end
