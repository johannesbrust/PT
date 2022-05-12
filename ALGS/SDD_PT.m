function [ M, out ] = SDD_PT( m, k, pars )
%SDD_PT "Shifted Diagonal Design" (for Pooling Tests)
%   
% Algorithm for generating the pooling matrix M, in which
% each of the 1:m*m individual samples 
%
%   * is in k tests
%   * each test contains exactly m different samples
%   * each pair of samples appears at most once.
%
% NOTE: The pool size m must be a prime number
%
% The matrix is k-1 disjunct and can exactly specify up to k-1 
% defectives in a pooling algorithm.
%
% Pooling test problem:
%
%   M * x = y, (goal to recover x given known M and y)
%
% where [t,n]=size(M), t<n, M is a binary matrix with tests (rows)
% and samples (sometimes called items) in its columns.
%
%
% INPUTS:
% m: Pool size
% k: Multiplicity (m*(k) tests required to exactly detect upto k-1 defectives)
% pars: Struct with optional parameters
%   pars.print: Flag to print infos
%
% OUTPUTS:
% M: Pooling matrix (binary, [m*k,m*m]=size(M))
% out: Struct with output values
%   out.time: Computational time
%--------------------------------------------------------------------------
% 03/30/21, J.B., Initial implementation
% 05/12/22, J.B., Preparation for release (and require prime pool size)

ctime = tic;

% Retrieve optional parameter values
if isfield(pars,'print')
    printInfo = pars.print;
else
    printInfo = 0;
end

if printInfo == 1

    fprintf('\n----------- SDD Algorithm ----------- \n'); % 41 chars    
    fprintf('Problem Size \n');
    fprintf('Pools (m):                     %i \n',m);
    fprintf('Multiplicity (k):              %i \n',k);
    fprintf('Expected tests (m*k):          %i \n',m*k);
    
end

% Check inputs
if (m+1)<k
    error('Input error: Need k<=m+1 \n');
end
% Check inputs
if isprime(m)==0
    error('Input error: Need prime number m \n');
end

% Construction of M
M = zeros(k*m,m*m,'int8');

%% Slopes (-1,-1/2,...,-1/(m-1))
for i=1:min((m-1),k) % Slopes
    slope = (i-1+m); % Determines shifts for indices
    for j=1:m % Pools
        cIdx = j; % Current index for sample
        for l=1:m % Pool elements
            % Check if sample falls into current line or needs to be 
            % modulated
%             if j+(l-1)*slope <= m*l
%                 cIdx = j+(l-1)*slope;
%             else
%                 cIdx = j+(l-2)*slope + mod(j+(l-1)*slope,m*l);
%             end
            if l>1
                if cIdx + slope <= m*l % cIdx + (l-1)*slope <= m*l                
                    cIdx = cIdx + slope; % falls in line
                %cIdx = cIdx + (l-1)*slope; % falls in line
                else
                    cIdx = m*(l-1) + mod(cIdx+slope,m*l); % modulated
                %cIdx = cIdx + mod(cIdx+(l-1)*slope,m*l); % modulated          
                end
            end
            M(j+m*(i-1),cIdx) = 1;
        end
    end
end

%% Horizontal and vertical cases (i.e., slopes 0, inf)
% Vertical (slope inf)
if m<=k
    for j=1:m
        cIdx=j;
        for l=1:m        
            if l>1
                if (j+1)==l
                    cIdx = cIdx + m+(m-1);
                else
                    cIdx = cIdx + (m-1);
                end            
            end
            M(j+m*(m-1),cIdx) = 1;
        end
    end
end

% Horizontal (slope 0)
if k==(m+1)
    for j=1:m        
        for l=1:m            
            %cIdx = (j-1)*m+l;
            M(j+m*m,(j-1)*m+l) = 1;            
        end
    end
end

out.time = toc(ctime);

if printInfo == 1
    
    fprintf('OUTPUT(S):############## \n');
    fprintf('Time (constr):        %i \n',out.time);    
    fprintf('\n');
end
