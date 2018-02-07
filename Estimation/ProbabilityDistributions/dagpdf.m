function y = dagpdf(x,p,a,b)
%DAGPDF Dagum probability density function
%   Y = DAGPDF(X,P,A,B) returns the probability density function of the
%   Dagum Distribution with shape parameters P and A and scale
%   parameter B, at the values in X.
%
%   The size of Y is the common sizes of the inputs. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        P, A , B > 0
%
%   Note: The Dagum Distribution is also known as the Burr III
%   Distribution, and is also a special case of the Generalized
%   Beta Prime Distribution.
%
%   See also DAGCDF, DAGINV, DAGSTAT, DAGFIT, 
%            DAGLIKE, DAGRND, DAGSF, DAGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 4
    error('dagpdf:TooFewInputs',...
        'Requires four input argument.');
end

try
    %Special case of the Generalized Beta Prime Distribution
    y=gbetaprpdf(x,p,1,a,b);
catch err
    error('dagpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Catch round off
y(y<0)=0;

end