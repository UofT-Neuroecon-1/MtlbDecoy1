function r = batesrnd(n,a,b,varargin)
%BATESRND Random arrays from the Bates Distribution
%   R = BATESRND(N,A,B) returns an array of random numbers chosen from
%   the Bates Distribution, of the mean of N independent and identically
%   distributed random variables uniformly distributed continuously
%   from A to B
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BATESRND(N,A,B,MM,NN,...) or R = BATESRND(N,A,B,[MM,NN,...]) 
%   returns an MM-by-NN-by-... array.
%
%   Type: Continuous, bounded
%   Restrictions:
%        A <= B
%        N > 1        (integer)
%
%   See also BATESPDF, BATESCDF, BATESINV, BATESSTAT, BATESFIT, BATESLIKE
%

%   Mike Sheppard
%   Last Modified 3-Jul-2011


if nargin < 3
    error('batesrnd:TooFewInputs','Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar
sz=cell2mat(varargin);

%Catch restrictions
a(a>b)=NaN; b(a>b)=NaN; n(n<0|n~=round(n))=NaN;

%Use direct method by definition of distribution
rnmat=unifrnd(a,b,[sz n]);
%Sum across last dimension
dim=length([sz n]);
r=sum(rnmat,dim);
%Divide by n for mean
r=r/n;

end