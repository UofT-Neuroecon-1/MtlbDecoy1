function h = fhaz(x,v1,v2)
%FHAZ   F hazard function.
%   H = FHAZ(X,V1,V2) returns the hazard function of the F distribution
%   with V1 and V2 degrees of freedom at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also FPDF, FCDF, FINV, FSTAT, FFIT, FLIKE, FRND, FSF
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 3
    error('fhaz:TooFewInputs',...
        'Requires three input arguments.');
end

try
    h = fpdf(x,v1,v2) ./ fsf(x,v1,v2);
catch
    error('fhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end