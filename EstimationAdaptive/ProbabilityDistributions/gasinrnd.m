function r = gasinrnd(alpha,a,b,varargin)
%GASINRND Random arrays from the Generalized Arcsine Distribution.
%   R = GASINRND(ALPHA,A,B) returns an array of random numbers chosen from
%   the Generalized Arcsine Distribution with shape parameter ALPHA
%   on the interval [A,B].
%
%   The size of R is the common size of the parameters if all are arrays.
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = GASINRND(A,B,M,N,...) or R = GASINRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, bounded
%   Restrictions:
%         A < B
%         0 < ALPHA < 1
%
%   See also GASINPDF, GASINCDF, GASININV, GASINSTAT, GASINFIT, GASINLIKE
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 3
    error('gasinrnd:TooFewInputs','Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (0<alpha & alpha<1);
alpha(~okparam)=NaN; a(~okparam)=NaN; b(~okparam)=NaN;

%Special case of betarnd, and will catch errors
r = betarnd(1-alpha,alpha,varargin{:});

%Transform to interval [a,b]
r = a + (b-a).*r;



end
