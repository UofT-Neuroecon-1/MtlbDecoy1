function s = gyulesf(x,a,b)
%GYULESF Generalized Yule–Simon survival function
%   S = GYULESF(X,A,B) returns the Generalized Yule-Simon survival 
%   function with shape parameter A and B, at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
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
%            GYULEFIT, GYULELIKE, GYULERND, GYULEHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012



if nargin < 3
    error('gyulesf:TooFewInputs',...
          'Requires at least three input arguments.'); 
end


try
    s = 1 - gyulecdf(x,a,b);
catch
    error('gyulesf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end