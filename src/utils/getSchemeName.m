function name = getSchemeName(flux_limiter,kappa)
% Decide scheme label based on flux limiter and kappa
if isscalar(flux_limiter)
    if flux_limiter ~= 18
        Name = ['Roe-MUSCL_kappa=' num2str(kappa)];
    else
        Name = 'Roe';
    end
elseif ~isscalar(flux_limiter) && any(flux_limiter==18)
    Name = ['Roe_v_Roe-MUSCL_kappa=' num2str(kappa)];
else
    Name = ['Roe-MUSCL_kappa=' num2str(kappa)];
end
end
