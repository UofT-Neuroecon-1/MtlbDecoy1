function r = maxwrnd(a,varargin)
%MAXWRND Random arrays from the Maxwell-Boltzmann Distribution
%   R = MAXWRND(A) returns an array of random numbers chosen the
%   Maxwell-Boltzmann Distribution with shape parameter A.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%   R = MAXWRND(A,MM,NN,...) or R = MAXWRND(A,[MM,NN,...]) 
%   returns an MM-by-NN-by-... array.

%   Mike Sheppard
%   Last Modified 24-Mar-2011


if nargin < 1
    error('maxwrnd:TooFewInputs','Requires at least one input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

sz=cell2mat(varargin);

rnmat=normrnd(0,a,[sz 3]);
%Sum the squares across last dimension
dim=length([sz 3]);
r=sum(rnmat.^2,dim);
%Square root
r=r.^(1/2);

end