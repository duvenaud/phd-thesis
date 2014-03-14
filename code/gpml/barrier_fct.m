function [f df] = barrier_fct(lh,covfunc,x,target)
% wrapper used during GP training: penalizes unreasonable signal-to-noise
% ratios
% 
% Marc Deisenroth, 2010-12-16

scale = 100;    % scaling factor
maxSNR = 1000;  % bound for SNR
order =  10;    % order of polynomial

[f df] = poly_barrier(lh, covfunc, x, target, maxSNR, order);

%%
function [f df] = poly_barrier(lh, covfunc, x, target, maxSNR, order)

% call gpr
[f df] = gpr(lh,covfunc,x,target);

D = length(lh)-2;

% signal-to-noise ratio
snr =  exp(lh(D+1)-lh(D+2));

% new objective function (penalizes too high SNRs)
f = f + (snr/maxSNR).^order;

% update corresponding gradients of the NLML
df(D+1) = df(D+1) + order*(snr/maxSNR).^(order-1)*snr/maxSNR;
df(D+2) = df(D+2) - order*(snr/maxSNR).^(order-1)*snr/maxSNR;

% % set max. length-scale ratio
% maxLSR = 1e5; % for the squared length-scales
% 
% % look at length-scales now
% [minLLS min_idx] = min(lh(1:D)); % smallest log-length-scale
% 
% for i = setdiff(1:D,min_idx) % don't change fct value with lsr = 1
%   
%   lsr = exp(2*lh(i)-2*minLLS); % length-scale-ratio; always > 1
%   
%   % new objective function (penalizes too high LSRs)
%   f = f + (lsr/maxLSR).^order; 
%   
%   % update corresponding gradients of the NLML
%   df(i) = df(i) + 2*order*(lsr/maxLSR).^(order-1)*lsr/maxLSR;
%   df(min_idx) = df(min_idx) - 2*order*(lsr/maxLSR).^(order-1)*lsr/maxLSR;
%   
% end




%% 
function [f df] = log_barrier(lh, covfunc, x, target, maxSNR, scale)
% for hard constraints

% call gpr
[f df] = gpr(lh,covfunc,x,target);

% new objective function
f = f - scale*log(maxSNR - exp(lh(end-1))/exp(lh(end)));

% derivatives
outer = -scale/(maxSNR - exp(lh(end-1))/exp(lh(end)));
inner1 = - exp(lh(end-1))/exp(lh(end));
inner2 =   exp(lh(end-1))/exp(lh(end));

df(end-1) = df(end-1) + outer*inner1;
df(end) =   df(end) + outer*inner2;
