function rnd = normratrnd(s1,s2,varargin)
%NORMRATRND Random arrays from the Normal Ratio Distribution
%   RND = NORMRATRND(S1,S2) returns an array of random numbers chosen from
%   the Normal Ratio Distribution, of the ratio of two Normal Distributions
%   with standard deviations S1 and S2.
%
%   The size of R is the common size of the parameters if all are arrays.
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   RND = NORMRATRND(S1,S2,M1,M2,...) or RND = NORMRATRND(S1,S2,[M1,M2,...])
%   returns an M1-by-M2-by-... array.
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('normratrnd:TooFewInputs',...
        'Requires at three two input argument.');
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN for out of range parameters.
s1(s1<0)=NaN; s2(s2<0)=NaN;

try
    scale=s1./s2;
    rnd = cauchyrnd(0,scale,varargin{:});
catch
    error('normratpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end