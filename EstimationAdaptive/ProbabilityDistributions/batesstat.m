function [m,v] = batesstat(n,a,b)
%BATESSTAT Mean and variance for the Bates Distribution
%   [M,V] = BATESSTAT(n,a,b) returns the mean and variance of the 
%   Bates Distribution, of the mean of n independent random variables
%   uniformly distributed from A to B.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        N > 1  (integer)
%
%   See also BATESPDF, BATESCDF, BATESINV, BATESFIT, 
%            BATESLIKE, BATESRND, BATESSF, BATESHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if (nargin<1)
    error('batesstat:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin==1
    a=0; b=1;
elseif nargin==2
    b=1;
end


[errorcode n a b] = distchck(3,n,a,b);

if errorcode > 0
    error('batesstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Solve for mean and variance
m=(a+b)/2;
v=((b-a).^2)./(12*n);

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>1 & n==round(n));
m(~okparam)=NaN;
v(~okparam)=NaN;

end