function [m,v]=moyalstat(u,s)
%MOYALSTAT Mean and variance of the Moyal Distribution
%   [M,V] = MOYALSTAT(U,S) returns the mean and variance of the Moyal
%   Distribution with location parameter U and scale parameter S.
%
%   Type: Continuous, unbounded
%   Restrictions:
%      U any real number
%      S>0
%
%   The size of the outputs is the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011



if nargin < 2
    error('moyalstat:TooFewInputs',...
        'Requires two input arguments.');
end

% Return NaN for out of range parameters.
s(s<0)=NaN;

%Constants
eulergamma=0.577215664901532860606512090082402431042159335939923598805;
log2=0.693147180559945309417232121458176568075500134360255254120;
pi2d2=4.9348022005446793094172454999380755676568497036203953132066;
    
try
    m=u+(s.*(eulergamma+log2));
    v=pi2d2.*(s.^2)+zeros(size(u)); % expand v's size to match u if necessary
catch
    error('moyalstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end