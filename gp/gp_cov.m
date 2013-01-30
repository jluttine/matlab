% GP_COV - This is a general documentation file for using covariance
%          functions.
%
% The covariance functions are used similarly. The covariance function is
% initialized by giving some function specific inputs, for instance,
%
%   COVFUNC = GP_COV_EXAMPLE(P1, P2, ...)
%
% where the number of inputs and their meaning depends on the covariance
% function. The input arguments can be, for instance, a (squared) distance
% matrix between the inputs, the dimensionality of the input space or
% another covariance function. Thus, be careful especially when changing
% covariance functions that you also change the inputs as needed! For
% instance, some covariance functions may use a distance matrix while
% some others a squared distance matrix.
%
% The returned covariance function COVFUNC can be used to obtain the
% covariance matrix by giving a parameter vector THETA:
%
%   K = COVFUNC(THETA)
%   [K, DK] = COVFUNC(THETA)
%
% where K is the covariance matrix and DK is a cell array of derivatives of
% K with respect to the parameters THETA. The length of THETA and the
% meaning of its elements depends on the covariance function. The parameters
% can be, for instance, characteristic length scale, magnitude or
% period. Again, be sure that you know which element of THETA corresponds to
% which parameter of the covariance function. Especially, if you modify your
% covariance function construction, also modify THETA appropriately!
%
% You can form complex covariance functions by using some functions to chain
% simple covariance functions. For instance, GP_COV_SCALE can be used to add
% a scaling parameter, GP_COV_SUM to sum two covariance functions and
% GP_COV_PRODUCT to multiply two covariance functions. Make sure how the
% additional parameters are added to the vector THETA and how the parameter
% vectors of multiple covariance functions are combined.
%
% For most covariance functions, you can fix some parameters:
%
%   COVFUNC = GP_COV_EXAMPLE(..., 'PARAMETER', VALUE)
%
% Then the parameter is not read from the vector THETA anymore. For
% instance, if you construct a covariance function as
%
%   COVFUNC = GP_COV_PERIODIC(SQRT(SQ_DIST(X1,X2)),'WAVELENGTH',10)
%
% then THETA should contain only one element, the smoothness parameter.
%
% You can obtain some dimensionality information by calling COVFUNC
% without any inputs:
%
%   [N_THETA, N1, N2] = COVFUNC()
%
%   N_THETA : Number of non-fixed parameters
%   N1      : Number of rows in the covariance matrix
%   N2      : Number of columns in the covariance matrix
%
% In order to obtain the covariance matrix of a covariance function with no
% parameters, you should call COVFUNC([]), not COVFUNC(). Note the
% difference!
%
% At least the following covariance functions are implemented:
%
%   GP_COV_SE       : Isotropic squared exponential
%   GP_COV_RQ       : Isotropic rational quadratic
%   GP_COV_PP       : Piecewise polynomial
%   GP_COV_DELTA    : Delta function, isotropic noise
%   GP_COV_SCALE    : Scales a covariance function
%   GP_COV_SUM      : Sums two covariance functions
%   GP_COV_PRODUCT  : Multiplies two covariance functions
%   GP_COV_JITTER   : Adds a small constant to the diagonal for numerical
%                     stability
%   GP_COV_TOEPLITZ : Creates a covariance matrix with Toeplitz structure
%   GP_COV_WRAP     : Makes a covariance function from a covariance matrix
%
% You can write your own covariance functions as long as they fulfill the
% specification described here. In some cases, it might be a good idea to
% write some wrapper functions for the covariance functions to make the
% usage more straightforward. For instance,
%
%   COVFUNC_SE = @(THETA,X1,X2) FEVAL(GP_COV_SE(SQ_DIST(X1,X2)),THETA)
%
% would be used as
%
%   [K, DK] = COVFUNC_SE(THETA, X1, X2)
%
% which some users may prefer. However, this kind of interface may reduce
% efficiency, as one may need to compute same things several times (squared
% distance matrix in the above example).
%
% Note: Despite the similar naming convention, GP_COV_PSEUDO is
% fundamentally different kind of covariance function and can not be used
% similarly to other covariance functions. See GP_COV_PSEUDO for more
% information.
%                  
% See also GP_LEARN, GP_LEARN_PSEUDO, GP_PREDICT, GP_PREDICT_PSEUDO.

% Last modified 2010-01-27
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function gp_cov()

help gp_cov