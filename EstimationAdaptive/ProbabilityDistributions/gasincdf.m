function p = gasincdf(x,alpha,a,b)
%GASINCDF Generalized arcsine cumulative distribution function.
%   P = GASINCDF(X,ALPHA,A,B) returns the cumulative distribution function
%   of the Generalized Arcsine Distribution with shape parameter ALPHA
%   on the interval [A,B] at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1 respectively.
%
%   Distribution: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        0 < ALPHA < 1
%
%   See also GASINPDF, GASININV, GASINSTAT, GASINFIT, 
%            GASINLIKE, GASINRND, GASINSF, GASINHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if (nargin~=2)&&(nargin~=4)
    error('gasincdf:TooFewInputs',...
        'Requires either two or four input arguments.');
end

if nargin == 2
    a = 0;
    b = 1;
elseif nargin==3
    b=1;
end

[errorcode x alpha a b] = distchck(4,x,alpha,a,b);

if errorcode > 0
    error('gasincdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (-Inf<a & a<Inf) & (-Inf<b & b<Inf) & (a<b) & (0<alpha & alpha<1);

%Use BETACDF to catch additional errors
p = betacdf((x-a)./(b-a), 1-alpha, alpha);
p(~okparam)=NaN;

%Catch round off
p(p<0)=0; p(p>1)=1;

end