function x = wgnscrnd(R,varargin)
%WGNSCRND Random arrays from the Wigner semicircle distribution
%   X = WGNSCRND(R) returns an array of random numbers from the Wigner
%   semicircle distribution with radius R
% 
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   X = WGNSCRND(R,M,N,P,Q,...) or X = WGNSCRND(R,[M,N,P,Q,...])
%   returns a M-by-N-by-P-by-Q-by... array.

%   Mike Sheppard
%   Last Modified 22-Apr-2011

if nargin < 1
    error('wgnscrnd:TooFewInputs','Requires at least one input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

%Use wgnscinv
x = wgnscinv(rand(varargin{:}),R);

end

