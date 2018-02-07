function x = trnorminv(p,mu,sigma,xmin,xmax)
%TRNORMINV Truncated Normal Distribution
%   X = TRNORMINV(P,MU,SIGMA,XMIN,XMAX) returns the cumulative distribution of
%   the Truncated Normal Distribution with (untruncated) Normal
%   Distribution parameters of MU and SIGMA, truncated between XMIN and
%   XMAX
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 5
    error('trnorminv:TooFewInputs',...
          'Requires five input arguments.'); 
end

[errorcode p mu sigma xmin xmax] = distchck(5,p,mu,sigma,xmin,xmax);

if errorcode > 0
    error('trnorminv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(p,'single') || isa(mu,'single') || isa(sigma,'single') || ...
        isa(xmin,'single') || isa(xmax,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


z2=(xmax-mu)./sigma; z3=(xmin-mu)./sigma;
term1=normcdf(z2)-normcdf(z3);
term2=normcdf(z3)+(p.*term1);
x=mu+sigma.*norminv(term2);

x(p<0|p>1)=NaN;
%Catch boundary
x(p==0)=xmin(p==0); x(p==1)=xmax(p==1);



end

