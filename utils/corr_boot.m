function boot_ci = corr_boot(x, y, niter_cor_perm, ci_bounds)

if nargin == 2
    niter_cor_perm = 1000;
    ci_bounds = [.025, .975];
elseif nargin == 3
    ci_bounds = [.025, .975];
end

n = length(x);
boot_ci = zeros(1, niter_cor_perm);
for ii = 1 : niter_cor_perm
    idx = randsample(n, n, true);  
    boot_ci(ii) = corr(x(idx), y(idx));
end
boot_ci = sort(boot_ci);
boot_ci = boot_ci([floor(niter_cor_perm * ci_bounds(1)), floor(niter_cor_perm * ci_bounds(2))]);



end