function r = asinrnd(a,b,varargin)
%ASINRND Random arrays from the Arcsine Distribution.
%   R = ASINRND(A,B) returns an array of random numbers chosen from an
%   Arcsine Distribution on the interval [A,B].  The size of R is the 
%   common size of MU and SIGMA if both are arrays.  If either
%   parameter is a scalar, the size of R is the size of the other
%   parameter.
%
%   R = ASINRND(A,B,M,N,...) or R = ASINRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%         A < B
%
%   See also ASINPDF, ASINCDF, ASININV, ASINSTAT, 
%            ASINFIT, ASINLIKE, ASINSF, ASINHAZ
%

%   ASINRND uses the inversion method.

%   Mike Sheppard
%   Last Modified 25-May-2012


if nargin < 2
    error('asinrnd:TooFewInputs','Requires at least two input arguments.'); 
end


if isempty(varargin), varargin={1}; end %Scalar

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b);
a(~okparam)=NaN; 
b(~okparam)=NaN;

%Use asininv to catch all errors
r = asininv(rand(varargin{:}),a,b);

end
