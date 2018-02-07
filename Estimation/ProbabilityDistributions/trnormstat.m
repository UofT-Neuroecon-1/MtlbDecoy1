function [m,v]=trnormstat(mu,sigma,xmin,xmax)
%TRNORMSTAT Mean and variance of the Truncated Normal Distribution
%   [M,V]=TRNORMSTAT(MU,SIGMA,XMIN,XMAX) returns the mean and variance of
%   the Truncated Normal Distribution with (untruncated) Normal
%   Distribution parameters of MU and SIGMA, truncated between XMIN and
%   XMAX
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.


%   Mike Sheppard
%   Last Modified 5-Jun-2011

if nargin < 1
    error('trnormstat:TooFewInputs',...
          'Requires at least one input argument.'); 
end

[errorcode mu sigma xmin xmax] = distchck(4,mu,sigma,xmin,xmax);

if errorcode > 0
    error('trnormstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(mu,'single') || isa(sigma,'single') || isa(xmin,'single') || isa(xmax,'single')
    m=zeros(size(mu),'single');
else
    m=zeros(size(mu));
end
v=m;

zmin=(xmin-mu)./sigma; zmax=(xmax-mu)./sigma;
pmin=normpdf(zmin); pmax=normpdf(zmax);
cmin=normcdf(zmin); cmax=normcdf(zmax);
den=cmax-cmin;
term1=(pmin-pmax)./den;

%Mean
m=mu+(term1.*sigma);

%Variance
term2=((zmin.*pmin)-(zmax.*pmax))./den;
v=(sigma.^2).*(1+term2-(term1.^2));


end