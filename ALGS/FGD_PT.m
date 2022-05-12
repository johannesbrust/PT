function [ M, out ] = FGD_PT( p, n, k, pars )
%FGD_PT "Finite Geometry Design" (for Pooling Tests)
%
% Supports affine and projective geometries
%   
% Algorithm for generating the pooling matrix M, in which
% each of the 1:p^(2n) individual samples 
%
%   * is in at most k tests
%   * each test contains at amost m different samples
%   * each pair of samples appears at most once.
%
% NOTE: The pool size is a prime power m = p^n,
%   where p is prime and n > 0. For finite field computations 
%   this method uses an external function to construct finite fields
%
% The matrix is k disjunct and can exactly specify up to k 
% defectives in a pooling algorithm.
%
% Pooling test problem:
%
%   M * x = y, (goal to recover x given known M and y)
%
% where [t,n]=size(M), t<n, M is a binary matrix with tests (rows)
% and samples (sometimes called items) in its columns.
%
% INPUTS:
% p: Prime number
% n: exponent 
% k: Multiplicity (p^n*(k+1) tests required to exactly detect upto k defectives)
% pars: Struct with optional parameters
%   pars.print: Flag to print infos
%   pars.geo = 1 (affine geometry)
%   pars.geo = 2 (projective geometry)
%
% OUTPUTS:
% M: Pooling matrix (sparse, [m*k,m*m]=size(M))
% out: Struct with output values
%   out.time: Computational time
%--------------------------------------------------------------------------
% 05/12/22, J.B., Preparation for release 

ctime = tic;

% Retrieve optional parameter values
if isfield(pars,'print')
    printInfo = pars.print;
else
    printInfo = 0;
end
if isfield(pars,'geo')
    geo = pars.geo;
else
    geo = 1;
end

if printInfo == 1

    m = p^n;
    fprintf('\n----------- FGD Algorithm ----------- \n'); % 41 chars    
    fprintf('Problem Size \n');
    fprintf('Pools (m):                     %i \n',m);
    fprintf('Disjunct (k):                  %i \n',k);
    fprintf('Expected tests:                %i \n',m*(k+1));
    
end

% Check inputs
if (p^n+1)<k
    error('Input error: Need k <= "pool size"+1 \n');
end
% Check inputs
if isprime(p)==0
    error('Input error: Need prime number p \n');
end

pnum = 0;
if(n==1)
    pnum=1;
end

% Invoke a "Galois field".
% Constructing of pools is done along lines in the field
gfp=gf(p,n);

% Variables to setup lines
% Initializing slopes
s = 0:(p^n-1);
ls = length(s);

% Defining quantities for sparse representation
k1  = k+1;

XV = repmat(s',ls,1);
OV = ones(length(XV),1);
BV = reshape(repmat(s,ls,1),ls*ls,1);
SV = (0:ls*ls-1)';

% Sparse indices for the affine geometry
XIDX = zeros(ls*ls*k1,1); 
YIDX = zeros(ls*ls*k1,1); 

% Sparse indices for the projective geometry
if k < ls
    psize = ls-k;
else
    psize = ls+1;
end
XIDX1 = zeros(ls*(ls+1)*k1+psize,1); 
YIDX1 = zeros(ls*(ls+1)*k1+psize,1);

IDXlv2 = 0:ls*ls-1;
IDXlv = 0:ls-1;

for i=0:k-1
    
    % Sparse computations. Defined on vectors

    AV = s(i+1).*OV;
    
    if pnum==0
        AXV = gfp.dmult(AV,XV);
        YV = gfp.add(AXV,BV);
    else
        AXV = mod(AV.*XV,p);
        YV = mod(AXV+BV,p);
    end
    
    % Affine plane
    XIDX(i*ls*ls + IDXlv2 + 1) = i*ls + BV + 1;
    YIDX(i*ls*ls + IDXlv2 + 1) = SV(sub2ind([ls,ls],XV+1,YV+1)) + 1;
    
    % Projective Geometry
    % Based on appending to an affine plane
    XIDX1(i*ls*ls + i*ls + IDXlv2 + 1) = XIDX(i*ls*ls + IDXlv2 + 1);
    YIDX1(i*ls*ls + i*ls + IDXlv2 + 1) = YIDX(i*ls*ls + IDXlv2 + 1);
    
    XIDX1(i*ls*ls + i*ls + ls*ls + IDXlv + 1) = i*ls + IDXlv + 1;
    YIDX1(i*ls*ls + i*ls + ls*ls + IDXlv + 1) = ls*ls + i  + 1;
    
end


% Singular layer (i.e., infty slopes)
XIDX(ls*ls*k + IDXlv2 + 1) = ls*k + BV + 1;
YIDX(ls*ls*k + IDXlv2 + 1) = SV(sub2ind([ls,ls],BV+1,XV+1)) + 1;

% Singular layer for projective geometry (i.e., infty slopes)
XIDX1(ls*ls*k + k*ls + IDXlv2 + 1) = XIDX(k*ls*ls + IDXlv2 + 1);
YIDX1(ls*ls*k + k*ls + IDXlv2 + 1) = YIDX(k*ls*ls + IDXlv2 + 1);

XIDX1(ls*ls*k + k*ls + ls*ls + IDXlv + 1) = k*ls + IDXlv + 1;
YIDX1(ls*ls*k + k*ls + ls*ls + IDXlv + 1) = ls*ls + k  + 1;
    

% Singelton "projective geometry" (this geometry is less regular,
% which is why the matrix is extended with an identity block
% to ensure identifying upto k positives
if k < ls
    XIDX1(ls*ls*k + k*ls + ls*ls + ls + (1:(ls - k))) = ls*(k) + ls + (1:(ls-k));
    YIDX1(ls*ls*k + k*ls + ls*ls + ls + (1:(ls - k))) = ls*ls + ((k+1+1):(ls+1));
    tsize = ls-k;
else
    XIDX1(ls*ls*k + k*ls + ls*ls + ls + (1:(ls + 1))) = ls*(k) + ls + 1;
    YIDX1(ls*ls*k + k*ls + ls*ls + ls + (1:(ls + 1))) = ls*ls + (1:(ls+ 1));
    tsize = 1;
end

% Sparse matrices
if geo == 1    
    % Affine geometry
    M = sparse(XIDX,YIDX,1,(k+1)*ls,ls*ls); 
else
    % Projective geometry
    M = sparse(XIDX1,YIDX1,1,(k+1)*ls+tsize,ls*ls+ls+1); 
end


out.time = toc(ctime);

if printInfo == 1
    
    fprintf('OUTPUT(S):############## \n');
    fprintf('Time (constr):        %i \n',out.time);    
    fprintf('\n');
end
