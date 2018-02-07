function [m,v] = solitstat(N)
%SOLITINV Soliton Distribution
%   K = SOLITINV(p,n) returns the probability density of the
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 1
    error('solitstat:TooFewInputs',...
          'Requires one input arguments.');
end


% Initialize P to NaN.
if isa(N,'single')
    m=zeros(size(N),'single');
else
    m=zeros(size(N));
end
v=m;

eulergamma=-psi(1);
pn=psi(N);
m=(1./N)+eulergamma+pn;
m2=(-1)+eulergamma+(1./N)+N+pn;
v=m2-(m.^2);

end