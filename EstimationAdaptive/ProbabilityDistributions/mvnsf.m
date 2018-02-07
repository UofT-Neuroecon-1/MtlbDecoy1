function [s,err] = mvnsf(varargin)
%MVNSF Multivariate normal survival function
%   S = MVNSF(X) returns the survivel function of the multivariate
%   normal distribution with zero mean and identity covariance matrix,
%   evaluated at each row of X.  Rows of the N-by-D matrix X correspond to
%   observations or points, and columns correspond to variables or
%   coordinates.  S is an N-by-1 vector.
%
%   S = MVNSF(X,MU,SIGMA) returns the survival function of the
%   multivariate normal distribution with mean MU and covariance SIGMA,
%   evaluated at each row of X.  MU is a 1-by-D vector, and SIGMA is a
%   D-by-D symmetric, positive definite matrix.  MU can also be a scalar
%   value, which MVNSF replicates to match the size of X.  If the
%   covariance matrix is diagonal, containing variances along the diagonal
%   and zero covariances off the diagonal, SIGMA may also be specified as a
%   1-by-D matrix containing just the diagonal. Pass in the empty
%   matrix for MU to use its default value when you want to only specify
%   SIGMA.
%
%   The multivariate normal survival function at X is defined as the
%   probability that a random vector V, distributed as multivariate normal,
%   will fall outside the semi-infinite rectangle with lower limits defined by
%   X, i.e., Pr{V(1)>X(1), V(2)>X(2), ... V(D)>X(D)}.
%
%   S = MVNSF(XL,XU,MU,SIGMA) returns the multivariate normal survival 
%   function evaluated over the rectangle (hyper-rectangle for D>2) 
%   with lower and upper limits defined by XL and XU, respectively.
%
%   [S,ERR] = MVNSF(...) returns an estimate of the error in S.  For
%   bivariate and trivariate distributions, MVNSF uses adaptive quadrature on
%   a transformation of the t density, based on methods developed by Drezner
%   and Wesolowsky, and by Genz, as described in the references.  The default
%   absolute error tolerance for these cases is 1e-8.  For four or more
%   dimensions, MVNSF uses a quasi-Monte Carlo integration algorithm based on
%   methods developed by Genz and Bretz, as described in the references.  The
%   default absolute error tolerance for these cases is 1e-4.  For
%   univariate distributions and when SIGMA is specified as a diagonal,
%   numerical integration is not used and ERR returns NaN.
%
%   [...] = MVNSF(...,OPTIONS) specifies control parameters for the numerical
%   integration, used to compute Y when D > 1 and is not specified as a 
%   diagonal.  This argument can be created by a call to STATSET.  Choices 
%   of STATSET parameters are:
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
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3];
%      [X1,X2] = meshgrid(linspace(-1,3,25)', linspace(-3,1,25)');
%      X = [X1(:) X2(:)];
%      p = mvnsf(X, mu, Sigma);
%      surf(X1,X2,reshape(p,25,25));
%
%   See also MVNPDF, MVNCDF, MVNINV, MVNSTAT, MVNFIT, MVNLIKE, 
%            MVNRND, MVNHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

[y,err] = mvncdf(varargin{:});
s=1-y;

end