function [m,v] = hmtstat(n)
%HMTSTAT Mean and variance for the Heads-Minus-Tails Distribution
%   [M,V] = HMTSTAT(n) returns the mean and variance of the
%   Heads-Minus-Tails Distribution given a fair coin is tossed
%   2N number of times.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Discrete, bounded, {0,...,N}
%   Restrictions:
%        N >= 1         (integer)
%
%   See also HMTPDF, HMTCDF, HMTINV, HMTFIT, HMTLIKE, HMTRND, 
%            HMTSF, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 11-May-2011


if nargin < 1, 
    error('hmtstat:TooFewInputs','Requires one input arguments.'); 
end


% Initialize the mean and variance to NaN.
if isa(n,'single')
   m = zeros(size(n),'single');
else
   m = zeros(size(n));
end
v = m;


m=exp(log(2)+gammaln(n+0.5)-(1/2)*log(pi)-gammaln(n));
v=(2.*n)-(m.^2);


end