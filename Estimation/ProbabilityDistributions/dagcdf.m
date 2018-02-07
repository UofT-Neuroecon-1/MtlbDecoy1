function y = dagcdf(x,p,a,b)
%DAGCDF Dagum cumulative distribution function
%   Y = DAGCDF(X,P,A,B) returns the cumulative distribution function of
%   the Dagum Distribution with shape parameters P and A and scale
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
%   See also DAGPDF, DAGINV, DAGSTAT, DAGFIT, 
%            DAGLIKE, DAGRND, DAGSF, DAGHAZ
%

%   Mike Sheppard
%   Last Modified 7-May-2011

if nargin ~= 4
    error('dagcdf:TooFewInputs',...
        'Requires four input arguments.');
end

try
    %Special case of the Generalized Beta Prime Distribution
    y=gbetaprcdf(x,p,1,a,b);
    %Catch round off
    y(y<0)=0; y(y>1)=1;
catch err
    error('dagcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
