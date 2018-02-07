function [s,err] = mvtsf(varargin)
%MVTSF Multivariate t survival function
%   S = MVTSF(X,C,DF) returns the survival function of the multivariate
%   t distribution with correlation parameters C and degrees of freedom DF,
%   evaluated at each row of X.  Rows of the N-by-D matrix X correspond to
%   observations or points, and columns correspond to variables or
%   coordinates.  Y is an N-by-1 vector.
%
%   C is a symmetric, positive definite, D-by-D correlation matrix.  DF is a
%   scalar, or a vector with N elements.
%
%   Note: MVTSF computes the SF for the standard multivariate Student's t,
%   centered at the origin, with no scale parameters.  If C is a covariance
%   matrix, i.e. DIAG(C) is not all ones, MVTCDF rescales C to transform it
%   to a correlation matrix.  MVTCDF does not rescale X.
%
%   The multivariate t survival function at X is defined as the
%   probability that a random vector T, distributed as multivariate t, will
%   fall outside the semi-infinite rectangle with lower limits defined by X,
%   i.e., Pr{T(1)>X(1), T(2)>X(2), ... T(D)>X(D)}.
%
%   S = MVTSF(XL,XU,C,DF) returns the multivariate t survival function
%   evaluated over the rectangle with lower and upper limits defined by XL and
%   XU, respectively.
%
%   [S,ERR] = MVTSF(...) returns an estimate of the error in S.  For
%   bivariate and trivariate distributions, MVTSF uses adaptive quadrature on
%   a transformation of the t density, based on methods developed by Genz, as
%   described in the references.  The default absolute error tolerance for
%   these cases is 1e-8.  For four or more dimensions, MVTSF uses a
%   quasi-Monte Carlo integration algorithm based on methods developed by Genz
%   and Bretz, as described in the references.  The default absolute error
%   tolerance for these cases is 1e-4.
%
%   [...] = MVTSF(...,OPTIONS) specifies control parameters for the numerical
%   integration used to compute S.  This argument can be created by a call to
%   STATSET.  Choices of STATSET parameters are:
%
%         'TolFun'      - Maximum absolute error tolerance.  Default is 1e-8
%                         when D < 4, or 1e-4 when D >= 4.
%         'MaxFunEvals' - Maximum number of integrand evaluations allowed when
%                         D >= 4.  Default is 1e7.  Ignored when D < 4.
%         'Display'     - Level of display output.  Choices are 'off' (the
%                         default), 'iter', and 'final'.  Ignored when D < 4.
%
%   Example:
%
%      C = [1 .4; .4 1]; df = 2;
%      [X1,X2] = meshgrid(linspace(-2,2,25)', linspace(-2,2,25)');
%      X = [X1(:) X2(:)];
%      p = mvtsf(X, C, df);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVNPDF, MVNCDF, MVNINV, MVNSTAT, MVNFIT, MVNLIKE, 
%            MVNRND, MVNHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011


[y,err] = mvtcdf(varargin{:});
s=1-y;

end