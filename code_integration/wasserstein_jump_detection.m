% ------------------------------------------------------------
% Load the saved struct
load AreaDistributions.mat    % gives variable areaStruct

% Helper: extract the JOO value from the folder name
extractJOO = @(name) sscanf(name,'JBB%*f_JOO%f_JBO%*f');

% Sort by JOO ascending
[~,idx] = sort(arrayfun(@(s) extractJOO(s.folder), areaStruct));
areaStruct = areaStruct(idx);

% Pre‑allocate distance array
n = numel(areaStruct);
W = zeros(n-1,1);

% Simple 1‑Wasserstein distance for 1‑D histograms
for k = 1:n-1
    p = areaStruct(k  ).hist;    % probability vectors
    q = areaStruct(k+1).hist;
    W(k) = sum( abs( cumsum(p) - cumsum(q) ) );  % L1 of CDF diff
end

% Plot
figure;
plot(0:n-2, W,'o-','LineWidth',1.5);
xlabel('interval index  (between JOO_k and JOO_{k+1})');
ylabel('1‑Wasserstein distance');
title('Consecutive Wasserstein distances along JOO line');

% Find the largest jump
[~,maxIdx] = max(W);
fprintf('Largest jump between index %d and %d  (JOO ≈ %.3f → %.3f)\n',...
        maxIdx-1, maxIdx, ...
        extractJOO(areaStruct(maxIdx  ).folder), ...
        extractJOO(areaStruct(maxIdx+1).folder));
