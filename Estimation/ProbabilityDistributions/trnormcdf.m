function p = trnormcdf(x,mu,sigma,xmin,xmax)
%TRNORMCDF Truncated Normal Distribution
%   P = TRNORMCDF(X,MU,SIGMA,XMIN,XMAX) returns the cumulative distribution of
%   the Truncated Normal Distribution with (untruncated) Normal
%   Distribution parameters of MU and SIGMA, truncated between XMIN and
%   XMAX
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 5
    error('trnormcdf:TooFewInputs',...
          'Requires five input arguments.'); 
end

[errorcode x mu sigma xmin xmax] = distchck(5,x,mu,sigma,xmin,xmax);

if errorcode > 0
    error('trnormcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single') || ...
        isa(xmin,'single') || isa(xmax,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


z1=(x-mu)./sigma; z2=(xmax-mu)./sigma; z3=(xmin-mu)./sigma;

%In terms of the Normal Distribution
p=(normcdf(z1)-normcdf(z3))./(normcdf(z2)-normcdf(z3));

%Catch out of bounds and round off
p(x<xmin)=0; p(x>xmax)=1; 

end

