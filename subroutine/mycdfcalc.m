function [yCDF,xCDF,n,emsg,eid] = mycdfcalc(x,xname)
%CDFCALC Calculate an empirical cumulative distribution function.
%   [YCDF,XCDF] = CDFCALC(X) calculates an empirical cumulative
%   distribution function (CDF) of the observations in the data sample
%   vector X. X may be a row or column vector, and represents a random
%   sample of observations from some underlying distribution.  On
%   return XCDF is the set of X values at which the CDF increases.
%   At XCDF(i), the function increases from YCDF(i) to YCDF(i+1).
%
%   [YCDF,XCDF,N] = CDFCALC(X) also returns N, the sample size.
%
%   [YCDF,XCDF,N,EMSG,EID] = CDFCALC(X) also returns an error message and
%   error id if X is not a vector or if it contains no values other than NaN.
%
%   See also CDFPLOT.

%   Copyright 1993-2004 The MathWorks, Inc. 


% Ensure the data is a VECTOR.
yCDF = [];
xCDF = [];
if (nargin < 2)
   if isempty(inputname(1))
      xname = getString(message('stats:cdfcalc:DefaultInputName'));
   else
      xname = inputname(1);
   end
end
n = 0;
if (min(size(x)) ~= 1)
    warning(message('stats:cdfcalc:VectorRequired', xname));
    emsg = getString(message('stats:cdfcalc:VectorRequired', xname));
    eid = 'VectorRequired';
    return
end

% Remove missing observations indicated by NaN's.
x = x(~isnan(x));
n = length(x);
if n == 0
    warning(message('stats:cdfcalc:NotEnoughData', xname));
    emsg = getString(message('stats:cdfcalc:NotEnoughData', xname));
    eid = 'NotEnoughData';
    return
end

% Sort observation data in ascending order.
x = sort(x(:));

%
% Compute cumulative sum such that the sample CDF is
% F(x) = (number of observations <= x) / (total number of observations).
% Note that the bin edges are padded with +/- infinity for auto-scaling of
% the x-axis.
%

% Get cumulative sums
yCDF = (1:n)' / n;

% Remove duplicates; only need final one with total count
notdup = ([diff(x(:)); 1] > 0);
xCDF = x(notdup);
yCDF = [0; yCDF(notdup)];
emsg = '';
eid = '';

xCDF = [-Inf; xCDF; Inf];
yCDF = [yCDF; 1];

