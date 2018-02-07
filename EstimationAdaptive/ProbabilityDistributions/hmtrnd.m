function r = hmtrnd(n,varargin)
%HMTRND Random arrays from the Heads-Minus-Tails distribution
%   R = HMTRND(N) returns an array of random numbers chosen from the
%   Heads-Minus-Tails distribution, given a fair coin is tossed 2N number
%   of times.

%   The size of R is the common size of the input arrays.  If any 
%   parameter is a scalar, the size of R is the size of the other parameters.
%
%   R = HMTRND(N,P,Q,...) or R = HMTRND(N,[P,Q,...]) returns
%   an P-by-Q-by-... array.
%

%   Mike Sheppard
%   Last Modified 13-May-2011



if nargin < 1
    error('hmtrnd:TooFewInputs','Requires at least one input arguments.');
end

if isempty(varargin), varargin={[1] [1]}; end %Scalar

sz=cell2mat(varargin);

rnmat=bernrnd(.5,[sz 2*n]);

%Heads-Minus-Tials across last dimension
dim=length([sz 2*n]);
r=abs(sum(rnmat==0,dim)-sum(rnmat==1,dim));

end
