function [f,c_len , ma_len, iota, res] = compute_f(objective_data, coilData, magnetic_axis_data, adjoint)    
    %% update all coil points and tangents on coils
    coilData.update();
    %% compute sigma and iota
    outres = magnetic_axis_data.determineSigmaAndIota();
    %% update QS magnetic field and gradient of magnetic field
    magnetic_axis_data.compute_QS_axis_data();
    %% update coil magnetic field and gradient of magnetic field
    magnetic_axis_data.compute_coil_axis_data(coilData, adjoint);

    c_len = coilData.get_clength();
    c_len_ind = zeros(1, coilData.C);
    for k = 1:coilData.C
       c_len_ind(k) = coilData.get_clength_ind(k); 
    end
    
    
    ma_len = magnetic_axis_data.get_alength();
    iota = magnetic_axis_data.iota;
    Bdiff = magnetic_axis_data.computeBdifference();
    gradBdiff = magnetic_axis_data.computeGradBdifference();
    res = max(outres);
    
    len_pen = 0;
    for k = 1:coilData.C
       len_pen = len_pen + 0.5*((c_len_ind(k) - objective_data.c_target)/objective_data.c_target)^2.; 
    end
    
    f = Bdiff + gradBdiff +...
      + len_pen ...
      + 0.5*((ma_len - objective_data.ma_len_target)/objective_data.ma_len_target)^2....
      + 0.5*((iota - objective_data.iota_target)/objective_data.iota_target)^2 ;
    
    objective_data.current_data = [f,c_len_ind, ma_len,iota, Bdiff, gradBdiff];
    objective_data.current_g = f;
    objective_data.current_c_len = c_len;
    objective_data.current_c_len_ind = c_len_ind;
    objective_data.current_ma_len = ma_len;
    objective_data.current_iota = iota;
end
