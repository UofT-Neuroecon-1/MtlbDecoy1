function rnd = normsumrnd(m,v,w,varargin)
%NORMSUMRND Random arrays from the Normal Sum Distribution
%   RND = NORMSUMRND(M,V,W) returns an array of random numbers chosen the
%   Normal Sum Distribution , the sum of Normal Distributions with Means M,
%   variances V, and weight W. The size of RND is the common size of the
%   input arguments if both are arrays. If any parameter is a scalar,
%   the size of RND is the common size of the other parameters.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   RND = NORMSUMRND(M,V,W,M1,M2,...) or RND = NORMSUMRND(M,V,W,[M1,M2,...])
%   returns an M1-by-M2-by-... array.
%

%   Mike Sheppard
%   Last Modified 18-Jun-2011


if nargin < 3
    error('normsumrnd:TooFewInputs',...
          'Requires at three two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

[errorcode, m,v,w] = distchck(3,m,v,w);

if errorcode > 0
    error('normsumpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

mu=dot(m,w);
vr=dot(v,w.^2);
sig=sqrt(vr);

rnd = normrnd(mu,sig,varargin{:});

end
