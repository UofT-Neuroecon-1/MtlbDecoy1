function p = batescdf(x,n,a,b)
%BATESCDF Bates cumulative distribution function
%   P = BATESCDF(x,n,a,b) returns the cumulative distribution function 
%   of the Bates Distribution of the mean of n independent random 
%   variables uniformly distributed from A to B.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        N > 1    (integer)
%
%   Note: Loss of accuracy occurs when N is greater than 25.
%
%   See also BATESPDF, BATESINV, BATESSTAT, BATESFIT, BATESLIKE, BATESRND, 
%            BATESSF, BATESHAZ
%

%   Mike Sheppard
%   Last Modified: 15-Dec-2011


if (nargin~=2)
    error('batescdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end

if nargin==2
    a=0; b=1;
elseif nargin==3
    b=1;
end


[errorcode x n a b] = distchck(4,x,n,a,b);

if errorcode > 0
    error('batescdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Transformation of Irwin-Hall Distribution
p = irwinhallcdf(x.*n,n,a,b);

% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (n>1 & n==round(n));
%okvar = (a <= x & x <= b);

p(~okparam)=NaN;
p(okparam & x<a)=0;
p(okparam & x>b)=1;

%Catch round off
p(p<0)=0; p(p>1)=1;

end