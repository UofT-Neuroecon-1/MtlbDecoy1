function rnd = pascalrnd(n,p,varargin)
%PASCALRND Random arrays from the Pascal Distribution
%   RND = PASCALRND(N,P) returns an array of random numbers chosen from a
%   Pascal Distribution with parameters N and P. The size of RND is the
%   common size of N and P if both are arrays. If either parameter is a
%   scalar, the size of RND is the size of the other parameter.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   RND = PASCALRND(N,P,M1,M2,...) or RND = NBINDRND(N,P,[M1,M2,...])
%   returns an M1-by-M2-by-... array.
%
%   PASCALRND uses a transformation of NBINRND, which uses either a sum of
%   geometrix random values, or a Poisson/gamma mixture.

%   Mike Sheppard
%   Last Modified 18-Jun-2011


if nargin < 2
    error('pascalrnd:TooFewInputs',...
          'Requires at least two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%PASCALPDF(X,N,P)=NBINPDF(X-N,N,P)
%PASCALPDF(X+N,N,P)=NBIN(X,N,P);
rnd = n+nbinrnd(n,p,varargin{:});


end
