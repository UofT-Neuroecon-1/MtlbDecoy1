function [m,v] = benistat(a,b,s)
%BENISTAT Mean and variance for the Benini distribution
%   [M,V] = BENISTAT(A,B,S) returns the mean and variance of the Benini
%   distribution with shape parameters A and B and scale parameter S.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   Note: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENICDF, BENIINV, BENIFIT, BENILIKE, BENIRND,
%            BENISF, BENIHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011



if nargin < 3
    error('benistat:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode a b s] = distchck(3,a,b,s);

if errorcode > 0
    error('benistat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<s & s<Inf);


%Mean calculate in terms
temp1=(-1+a)./(2.*sqrt(b));
temp2=(-2+a)./(2.*sqrt(b));
m=s+((exp(temp1.^2).*erfc(temp1).*sqrt(pi).*s)./(2*sqrt(b)));

%variance compute in terms
term1=4.*exp(temp2.^2).*sqrt(b).*erfc(temp2);
term2=4.*exp(temp1.^2).*sqrt(b).*erfc(temp1);
term3=exp(2*(temp1.^2)).*sqrt(pi).*(erfc(temp1)).^2;
v=sqrt(pi).*s.^2.*(term1-term2-term3)./(4*b);

%Catch out of bounds
m(~okparam)=NaN;
v(~okparam)=NaN;

end