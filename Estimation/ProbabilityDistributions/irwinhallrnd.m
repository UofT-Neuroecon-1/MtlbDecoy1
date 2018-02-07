function r = irwinhallrnd(n,a,b,varargin)
%IRWINHALLRND Random arrays from the Irwin-Hall Distribution
%   R = IRWINHALLRND(N,A,B) returns an array of random numbers chosen from
%   the Irwin-Hall distribution as the sum of N independent and identically
%   distributed random variables uniformly distributed continuously
%   from A to B
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = IRWINHALLRND(N,A,B,MM,NN,...) or R = IRWINHALLRND(N,A,B,[MM,NN,...]) 
%   returns an MM-by-NN-by-... array.
%
%   Type: Continuous, bounded, (A*N,B*N)
%   Restrictions:
%        A < B
%        N >= 1        (integer)
%
%   Note: The Irwin-Hall Distribution is also known as the
%   Uniform Sum Distribution. 
%
%   See also IRWINHALLPDF, IRWINHALLCDF, IRWINHALLINV, IRWINHALLSTAT, 
%            IRWINHALLFIT, IRWINHALLLIKE, IRWINHALLSF, IRWINHALLHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin < 3
    error('irwinhallrnd:TooFewInputs','Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

sz=cell2mat(varargin);

if isscalar(sz), sz=[sz sz]; end %Square matrix

%By definition
rnmat=unifrnd(a,b,[sz n]);
%Sum across last dimension
dim=length([sz n]);
r=sum(rnmat,dim);

end