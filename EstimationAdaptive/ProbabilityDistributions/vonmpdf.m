function y = vonmpdf(x,u,k)
%VONMPDF Von-Mises probability density function
%   Y = VONMPDF(X,U,K) returns the Von-Mises probability density function
%   with mean U and concentration K.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 30-May-2011


if nargin < 3
    error('vonmpdf:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode x u k] = distchck(3,x,u,k);

if errorcode > 0
    error('vonmpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


num=exp(k.*cos(x-u));
den=2*pi*besseli(0,k);
y=num./den;
%If u=0 it is equivalent to the uniform distribution
y(u==0)=1/(2*pi);
%Valid region only
y(x<u-pi | x>u+pi)=0;


end