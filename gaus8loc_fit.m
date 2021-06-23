function [R2 s m sem betas modelFun] = gaus8loc_fit(inp_rate,colors)
% Fits 1-D Gaussian curve to 8 spatial locations
% Input rates should be rotated so that:
% the best rate corresponds to location 5
% rate 9 is the same as rate 1
% rate 10 corresponds to the fovea
%
% Calls fittinggaus function
% Author:           Christos Constantinidis, Ph.D.
% Created:          11-Mar-2010

loc=1:9;
xgrid=1:0.1:9;
if min(size(inp_rate)) > 1
    [s m sem nt] = allmeans_xq([1:size(inp_rate,1)],inp_rate);
end
% rate(1:9)=inp_rate(1:9);
rate(1:9) = m(1:9);

modelFun =  @(p,x) p(1)+p(2) .* exp(-0.5 *((x - p(3))/p(4)) .^2);
beta_init = [min(rate) max(rate)-min(rate) 4 2];

% beta_init = [min(rate) max(rate)-min(rate) 5.03 0.8];
% beta_init = [min(rate) max(rate)-min(rate) 4 1];
% beta_init = [5.4 2.2 5.03 0.2]; % for pre vlpfc spatial task

betas = nlinfit(loc,rate,modelFun,beta_init);
rate_fit = modelFun(betas,loc);
correl_mat = corrcoef(rate,rate_fit);
R2 = correl_mat(2)*correl_mat(2);
betas = betas';
% rate_plot = modelFun(betas,xgrid);

% figure; 
% hold on;
% errorbar(s,m,sem,'marker','o','linestyle','none','color',colors)%(:,2)-m)
% plot(xgrid,rate_plot,colors);
% xlim([0.5 10.5]);

% 2d contourf plot
%  figure;
%  hold on
% rates = [m(6) m(7) m(8);m(5) m(9) m(1); m(4) m(3) m(2)];
% [X,Y]=meshgrid(-10:10:10, -10:10:10);
% [XI YI]=meshgrid(-10:1:10, -10:1:10);
% rateint=interp2(X,Y,rates,XI,YI);
% contourf(XI,YI,rateint,40)
% shading flat

%============================= end allmeans ===============================

function [s, m, sem, nt] = allmeans_xq(st, r)

% Given the stimulus (st) and the responses (r) at each trial -- trials
% correspond to rows in r -- this function returns a vector of stimuli
% without repeats and means, SEMs and numbers of trials. 

 % 
 % first make sure it's one neuron only and de-NaN, just to be safe
 %
 st = colvec(st);  
%  r  = colvec(r); % xq
try
 ii = find(isnan(st)==1 | isnan(r)==1);
 st(ii) = [];
 r(ii)  = [];
end
 %
 % put a minimum number of trials as a condition
 %
 [s nt] = nunique(st);
 %
 % compute means and sems at each stimulus value
 %
 Ns = length(s);
%  m  = zeros(Ns,1);
%  sem= zeros(Ns,1);
for j=1:Ns
    %      ii = find(r(j,:)~=nan);
    ii = find(~isnan(r(j,:)));
    m(j)  = nanmean(r(j,:));
    sem(j) = std(r(j,ii))./sqrt(length(ii));     % this can be zero
    nt1(j) = length(ii);
end
 nt = min(nt1);
 m = colvec(m);
 sem = colvec(sem);
%  sem = sem./sqrt(nt);

%============================= end allmeans ===============================