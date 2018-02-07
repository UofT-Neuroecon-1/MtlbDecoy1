function y = batespdf(x,n,a,b)
%BATESPDF Bates probability density function
%   Y = BATESPDF(x,n,a,b) returns the probability density function of the 
%   Bates Distribution representing the mean of N independent random
%   variables uniformly distributed from A to B.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        N > 1     (integer)
%
%   Note: Loss of accuracy occurs when N is greater than 25.
%
%   See also BATESCDF, BATESINV, BATESSTAT, BATESFIT, 
%            BATESLIKE, BATESRND, BATESSF, BATESHAZ
%

%   Mike Sheppard
%   Last Modified: 13-May-2012


if (nargin<2)
    error('batespdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end

if nargin==2
    a=0; b=1;
elseif nargin==3
    b=1;
end


[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('batespdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of Irwin-Hall Distribution
y = n.*irwinhallpdf(x.*n,n,a,b);

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>1 & n==round(n));
okvar = (a <= x & x <= b);

y(~okparam)=NaN;
y(okparam & ~okvar)=0;

%Catch round off
y(y<0)=0;

end