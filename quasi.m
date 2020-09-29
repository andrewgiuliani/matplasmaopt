function [f,df, res] = quasi(x, coilData, magnetic_axis_data, objective_data, idxs)

    magnetic_axis_data.eta_bar = x(idxs.eta_bar_idx) ;
    magnetic_axis_data.cR = x(idxs.cR_idx)';
    magnetic_axis_data.sZ = x(idxs.sZ_idx)';
    coilData.I = x(idxs.I_idx)' / ( 4 * pi * 10^(-7) );
    coilData.coil_coeffs(:) = x(idxs.cc_idx)';
    
    [f, c_len , ma_len, iota, res] = compute_f(objective_data, coilData, magnetic_axis_data, true);
    df = df_adjoint(objective_data, coilData, magnetic_axis_data);
    

end

