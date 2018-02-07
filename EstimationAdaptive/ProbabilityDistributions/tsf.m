function s = tsf(x,v)
%TSF    Student's T survival function
%   S = TSF(X,V) computes the survival function of the Student's T 
%   distribution with V degrees of freedom, at the values in X.
%
%   The size of S is the common size of X and V. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   See also TPDF, TCDF, TINV, TSTAT, TFIT, TLIKE, TRND, THAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


if nargin < 2,
    error(message('tsf:TooFewInputs'));
end


try
    s = 1 - tcdf(x,v);
catch
    error('tsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end