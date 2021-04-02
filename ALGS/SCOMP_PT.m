function [ x, out ] = SCOMP_PT( M, y, pars )
% Sequential COMP (SCOMP)
%   
% Algorithm for pooling tests (e.g., mixed COVID samples) with
% problem representation:
%
%   M * x = y, (goal to recover x given known M and y)
%
% where [t,n]=size(M), t<n, M is a binary matrix with tests (rows)
% and patients (sometimes called items) in its columns.
%
% Procedure:
% This method first applies the same steps as COMP, and subsequently
% searches for definite defectives (DD) (after removing any items that 
% are known to be non-defective).
% Subsequently, a further loop checks and improves on remaining
% "unexplained" tests.
%
% INPUTS:
% M: Matrix (t x n)
% y: Measured outcomes (t x 1)
% pars: Struct with optional parameters
%   pars.print: Flag to print infos
%   pars.d: Number of defectives
%
% OUTPUTS:
% x: Computed solution
% out: Struct with output values
%   out.time: Computational time
%   out.d_final: Computed defectives
%   out.uex: Unexplained tests
%   out.nCds: Number of trial (items) candidates
%   out.cds: Indices of trial (items) candidates 
%
%
% Note: This method assumes that inputs M and y are logical.
%--------------------------------------------------------------------------
% 03/18/21, J.B., Initial implementation

ctime = tic;

% Problem data
[t,n] = size(M);

% Retrieve optional parameter values
if isfield(pars,'print')
    printInfo = pars.print;
else
    printInfo = 0;
end
if isfield(pars,'d')
    d = pars.d;
else
    d = -1;
end

if printInfo == 1

    fprintf('\n----------- SCOMP Algorithm ----------- \n'); % 41 chars    
    fprintf('Problem Size \n');
    fprintf('Tests:                 %i \n',t);
    fprintf('Patients:              %i \n',n);
    fprintf('Input Expct. Positive: %i \n',d);
    
end

x = zeros(n,1);
% Algorithm (column-wise)
% COMP
for i=1:n
    if sum(M(:,i)) > 0
        x(i) = (sum(M(:,i)) == sum(1 == y(M(:,i))));
    end
    %x(i) = (sum(M(:,i)) == sum(y(M(:,i))));
end

% Algorithm
% DD
M1 = M(y>0,x>0);
x(x>0) = sum( M1(sum(M1,2) == 1,:), 1 ) > 0;
%x(x>0) = (1 == sum(M(y>0,x>0),2));

% Algorithm SCOMP
d_comp = sum(x);
yp = y>0;
xp = x>0;
iyp = find(yp);
M1 = M(yp,xp);
% Unexplained tests
yu = (sum(M1,2)-1) < 0;
iter = 1;

% (Additional) candidates checked but not declared positive
cds = zeros(n-d_comp,1,'int8');
nCds = 0;

while ( sum(yu) > 0 ) && ( iter < (n-d_comp+1) )

    %[mv,mx] = max(sum( M(iyp(yu),:), 1 ));
    [mv,mx] = sort(sum( M(iyp(yu),:), 1 ),'descend');
    
    % Check possible candidates
    for i=1:length(mx)
        
        if mv(i) > 0 
            
            x(mx(i)) = 1;
            
            % Check whether previous decisions have changed when new
            % item is added. In this case candidate item is not stored            
            if sum(sum(M(y==0,x>0),2) > 0) > 0
                
                x(mx(i)) = 0;
                
                % Store trial candidates
                if sum(cds == mx(i)) == 0
                    
                    nCds = nCds+1;
                    cds(nCds) = mx(i);
                    
                end
            
            else
                
                % Break out of loop if item is added during
                % current iteration
                break;
                
            end
            
        end
        
    end
    
    % Updates
    xp = x>0;
    M1 = M(yp,xp);    
    yu = (sum(M1,2)-1) < 0;
    iter = iter+1;

end    

out.d_final = sum(x);
out.time = toc(ctime);
out.uex = sum(yu);
out.nCds = nCds;
out.cds = cds;

if printInfo == 1
    
    fprintf('OUTPUTS:################ \n');
    fprintf('Time (search):         %i \n',out.time);
    fprintf('Positive items:        %i \n',out.d_final);
    fprintf('No. unexpl. tests:     %i \n',out.uex);
    fprintf('Indices unexpl. tests: %i \n',iyp(yu));
    fprintf('No. tried cands:       %i \n',out.nCds);
    fprintf('Indices tried cands:   %i \n',out.cds(1:nCds));
    fprintf('######################## \n');
    fprintf('Identified Indices:    %i \n',find(x));
    fprintf('\n');
    
end
