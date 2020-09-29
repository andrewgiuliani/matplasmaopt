function [grad_g] = df_adjoint(objective_data, coilData, magnetic_axis_data)
    
    % assumes compute f has been executed
    magnetic_axis_data.update();
    
    % compute f_ca and f_s
    f_ca = magnetic_axis_data.compute_fca() ;
    f_s  = magnetic_axis_data.compute_fs() ;
    f_etabar = magnetic_axis_data.compute_fetabar();
    
    % compute g_ca
    g_ca =       (objective_data.current_ma_len-objective_data.ma_len_target) / (objective_data.ma_len_target^2) * magnetic_axis_data.compute_dalength_dca();
    g_ca = g_ca + magnetic_axis_data.compute_dBdiff_dca();
    g_ca = g_ca + magnetic_axis_data.compute_dGradBdiff_dca();

    
    % compute g_s
    g_s =       (objective_data.current_iota-objective_data.iota_target) / (objective_data.iota_target^2) * magnetic_axis_data.compute_iota_ds();
    g_s = g_s + magnetic_axis_data.compute_gradBdiff_ds();
    
    % compute g_cc
    % penalize each coil length separately
    indcoil_len = zeros(1, numel(coilData.coil_coeffs));
    vlen = numel(coilData.coil_coeffs)/coilData.C;
    
    for j = 1:numel(coilData.coil_coeffs)
        [k, kk] = ind2sub(size(coilData.coil_coeffs),j) ;
        i = ceil(kk/ 2/3);
        coord = mod(ceil(kk/2)-1,3)+1; %x,y,z
        
        indcoil_len(1,(i-1)*vlen+1:i*vlen) = (objective_data.current_c_len_ind(i)-objective_data.c_target) / (objective_data.c_target^2) * ones(1,vlen);
        if coord == 1 || coord == 2
            continue
        end
    end
            
                
    dg_dcc =  indcoil_len.* coilData.compute_dlen_dc();

    dg_dcc = dg_dcc + magnetic_axis_data.compute_dBdiff_dcc();
    dg_dcc = dg_dcc + magnetic_axis_data.compute_dGradBdiff_dcc();
    
 
    % g_etabar
    g_etabar = magnetic_axis_data.compute_dGradBdiff_detabar();
    
    % compute g_I
    dg_dI =         magnetic_axis_data.compute_dBdiff_dI();
    dg_dI = dg_dI + magnetic_axis_data.compute_dGradBdiff_dI();
    
    lambda = f_s'\g_s';
    
    dg_dca = g_ca - lambda'* f_ca;
    dg_etabar = g_etabar - lambda'* f_etabar;
    
    grad_g = [dg_etabar dg_dca dg_dI dg_dcc];
    
    objective_data.current_grad = grad_g;
    objective_data.current_gradnorm = max(abs(grad_g));
end


