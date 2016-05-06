% This script tests the minimization of non-convex quadratic program
%
% Written by Amit Kumar K C, PhD at Wollman Lab (UCSD)
% Date: 5th May 2015

n=100;
maxiter=200;
params=[1,n,1e-4, 100, 0, 0, 1, 1, 1, 1, 0];
solutions=zeros(n,maxiter,2);
fvals=zeros(maxiter,2);
for iter=1:maxiter

    c=randn(1);
    v=randn(n,1);

    %MATLAB's fmincon solution

    [y1,f1]=fmincon(@(x) (0.5*c*(x'*x)-v'*x), ones(n,1)/n, [],[],ones(1,n),1,zeros(n,1),[],[],optimset('Display','off','Algorithm','active-set'));


    y=mexSimplex(v,c,params);
    f=(0.5*c*(y'*y)-v'*y);

    solutions(:,iter,1)=y1;
    solutions(:,iter,2)=y;

    fvals(iter,1)=f1;
    fvals(iter,2)=f;

    prev_prog=curr_prog;
    curr_prog=floor(iter*100/maxiter);
    if mod(curr_prog, 5)==0 && curr_prog~=prev_prog
        fprintf('  %d  ', curr_prog);
    end
end
fprintf('\n');

%%
figure(1); plot(fvals(:,1), fvals(:,2),'r.'); title('functional values'); xlabel('matlab'); ylabel('mexSimplex');
figure(2); hold on;
for iter=1:maxiter
    plot(solutions(:,iter,1), solutions(:,iter,2),'-', 'color',rand(3,1));
end
hold off;
drawnow;
