function p = warcdf(x,a,b)
%WARCDF Waring cumulative distribution function
%   P = WARCDF(X,A,B) returns the Waring cumulative distribution
%   function with shape parameter A, at the values in X.
%
%   Note: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 25-May-2011



if nargin < 3
    error('gyulepdf:TooFewInputs',...
          'Requires at least three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('gyulepdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Uses the Pochhammer symbol, which is Poch[n,x]=Gamma[n+x]/Gamma[n]. Use
%gammaln for accuracy
num=gammaln(a+b)+gammaln(1+b+floor(x));
den=gammaln(b)+gammaln(1+a+b+floor(x));
p=1-exp(num-den);



end