function h = bortanhaz(x,alpha,n)
%BORTANHAZ Borel-Tanner hazard function
%   H = BORTANHAZ(X,A,N) returns the hazard function of the 
%   Borel-Tanner Distribution with shape parameters ALPHA and N,
%   at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        1 <= N         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANPDF, BORTANCDF, BORTANINV, BORTANSTAT, 
%            BORTANFIT, BORTANLIKE, BORTANRND, BORTANSF
%


%   Mike Sheppard
%   Last Modified 21-Dec-2011



if nargin<3
    error('bortanhaz:TooFewInputs','Requires three input argument.');
end


try
    yt = bortanpdf(x,alpha,n);
    st = bortansf(x,alpha,n);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('bortanhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end