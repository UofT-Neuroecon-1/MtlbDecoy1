function chihat = argusfit(x,c)
%ARGUSFIT Parameter estimate for ARGUS data
%   CHIHAT = ARGUSFIT(X,C) returns the maximum likelihood estimates of the
%   parameter CHI of the ARGUS distribution given the data in X and value
%   for C.
%
%   Default value for C is 1.
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSSTAT, 
%            ARGUSLIKE, ARGUSRND, ARGUSSF, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified: 2-Dec-2011


if nargin < 1
    error('argusfit:TooFewInputs',...
        'Requires at least one input arguments.');
end

if nargin == 1
    c=1;
end

if ~isscalar(c)||(c<0)||~isfinite(c)
    error('argusfit:TooFewInputs',...
        'Second input, C, must be a scalar valued non-negative real number.');
end


%MLE is the solution to the non-linear equation
%1 - (3/chi^2) + chi*phi(chi) / psi(chi) = (1/n)*sum(x_i^2/c^2)  [1]
%Where phi(x) is normpdf(x)
%and psi(x)=normcdf(x) - x*normpdf(x) - 1/2;

vals=(x(:)./c(:)).^2;
RHS=mean(vals); %RHS of [1]

if RHS<.4
    error('argusfit:SolveConditions',...
        'Unique solution does not exist.');
elseif RHS>1
    error('argusfit:SolveConditions',...
        'Solution does not exist for value of C given');
end

% Use FZERO to estimate solution
% f(chi_0) = LHS(chi_0) - RHS = 0
% In order for it to converge properly, use known values to interpolate
% good starting point.
lookup=[.4 9.933849125615400e-006;
        .5 1.676095720944136;
        .6 2.387951648352959;
        .7 3.040660772465421;
        .8 3.855439870439049;
        .9  5.477188957058711;
        .99  17.320508075688721];
chi_0 = interp1(lookup(:,1),lookup(:,2),RHS);    
if isnan(chi_0), chi_0=20; end %beyond RHS=.99
%Try fzero first
[chihat,fval1,exitflag1] = fzero(@(x)(LHS(x)-RHS),chi_0);
if exitflag1~=1
    %Try fminsearch
    [chihat,fval2,exitflag2]=fminsearch(@(x) (LHS(x)-RHS)^2,chi_0);
    if exitflag2~=1
        error('argusfit:SolveConditions',...
        'Solution not found');
    end
end

%LHS is an even function, as last error check in case negative solution was
%found, take absolute value
chihat=abs(chihat);

end



function y=LHS(chi)
psi=normcdf(chi)-chi.*normpdf(chi) - (1/2);
term1=1;
term2=-3./(chi.^2);
term3=chi.*normpdf(chi)./psi;
y = term1+term2+term3;
end


