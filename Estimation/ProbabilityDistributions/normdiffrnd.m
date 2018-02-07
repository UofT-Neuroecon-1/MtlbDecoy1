function rnd = normdiffrnd(m,v,varargin)
%NORMDIFFRND Random arrays from the Normal Difference Distribution
%   RND = NORMDIFFRND(M,V) returns an array of random numbers chosen from
%   the Normal Difference Distribution, the difference of two Normal
%   Distributions with means M, and variances V. 
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   RND = NORMDIFFRND(M,V,M1,M2,...) or RND = NORMDIFFRND(M,V,[M1,M2,...])
%   returns an M1-by-M2-by-... array.
%

%   Mike Sheppard
%   Last Modified 19-Jun-2011


if nargin < 2
    error('normdiffrnd:TooFewInputs',...
          'Requires at two two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Check sizes
if ~isvector(m) || ~isvector(v) || numel(m)~=2 || numel(v)~=2
    error('normdiffrnd:TooFewInputs',...
          'Mean and variances must be vectors of length 2.');
end

m=m(:)'; v=v(:)';
w=[1 -1]; %Use Normal Sum with weights of +1 and -1
rnd=normsumrnd(m,v,w,varargin{:});

end