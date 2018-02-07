function h = gyulehaz(x,a,b)
%GYULEHAZ Generalized Yule-Simon hazard function
%   H = GYULEHAZ(X,A,B) returns the hazard function of the Generalized
%   Yule-Simon Distribution with shape parameter A and B, at the values
%   in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Discrete, semi-bounded, {1,...,Inf}
%   Restrictions:
%      A > 0
%      0 <= B < 1
%
%   NOTE: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1. The two-parameter generalization of the
%   original Yule distribution replaces the beta function with an
%   incomplete beta function. 
%
%   See also GYULEPDF, GYULECDF, GYULEINV, GYULESTAT, 
%            GYULEFIT, GYULELIKE, GYULERND, GYULESF
%


%   Mike Sheppard
%   Last Modified 13-May-2012



if nargin<3
    error('gyulehaz:TooFewInputs','Requires three input argument.');
end


try
    yt = gyulepdf(x,a,b);
    st = gyulesf(x,a,b);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('gyulehaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end