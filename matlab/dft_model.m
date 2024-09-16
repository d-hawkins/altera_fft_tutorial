% -----------------------------------------------------------------------------
% dft_model.m
%
% 9/1/2024 D. W. Hawkins (dwh@caltech.edu)
%
% DFT model using double-precision calculations.
%
% -----------------------------------------------------------------------------
function X = dft_model(x)

	% Length
	N = length(x);

	% Calculate the frequency response
	X = zeros(N,1);
	n = [0:N-1]';
	for k = 0:N-1,
		X(k+1) = sum(x.*exp(-j*2*pi*n*k/N));
	end
end
