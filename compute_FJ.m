function [F,J] = compute_FJ(state, eta_bar, curvature, torsion, Dtilde, I2_over_B0, abs_G_over_B0)
    sigma = state(1:end-1);
    iota = state(end);
    
    F = zeros(size(state,1),1);
    F(1:end-1) =  Dtilde * sigma ...
                 + iota  * ( (eta_bar*eta_bar*eta_bar*eta_bar) ./ (curvature .* curvature .* curvature .* curvature) + 1 + sigma .* sigma) ...
                 - 2 * (eta_bar*eta_bar) ./ (curvature .* curvature) .* (I2_over_B0 - torsion) * abs_G_over_B0;
    F(end) = sigma(1) - 0.;
    
    
    J = zeros(size(state,1));
    J(1:end-1,end) =  (eta_bar*eta_bar*eta_bar*eta_bar) ./ (curvature .* curvature .* curvature .*curvature) + 1 + sigma .* sigma;
    J(1:end-1,1:end-1) = Dtilde + diag(iota * 2 * sigma);
    J(end, 1) = 1;
    
end

