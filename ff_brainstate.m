function [ff] = ff_brainstate(ff,brainstate)
% ff(:,2) is time stamps
% ff(:,1) is ff signal
idx = ceil((ff(:,2)+1e-13)/5);
idx2 = find(idx<length(brainstate)+1);
for i = 1:length(idx2)
    ff(i,3) = brainstate(idx(idx2(i)));
end
