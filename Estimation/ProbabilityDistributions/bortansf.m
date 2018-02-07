function s = bortansf(x,alpha,n)
%BORTANSF Borel-Tanner survival function
%   S = BORTANSF(X,A,N) returns the survival function of the 
%   Borel-Tanner Distribution with shape parameters ALPHA and N,
%   at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        1 <= N         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANPDF, BORTANCDF, BORTANINV, BORTANSTAT, 
%            BORTANFIT, BORTANLIKE, BORTANRND, BORTANHAZ
%


%   Mike Sheppard
%   Last Modified 21-Dec-2011



if nargin<3
    error('bortansf:TooFewInputs','Requires three input argument.');
end


try
    s = 1 - bortancdf(x,alpha,n);
catch
    error('bortansf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end