    gamma = 4 * pi * 10^(-7) ;
    c_target = 2*pi * 0.7;
    ma_len_target = 2 * pi;
    iota_target = 0.103;


    objective_data = objectiveDataType(c_target,ma_len_target,iota_target);
    %% Read in the initial coils
    coilData = coilDataType('./input/coils.dat', './input/I.dat', 2, 1, 10);
    magnetic_axis_data = magneticAxisType(coilData, './input/cR.dat',  './input/sZ.dat',   10); 
    

    eta_bar_idx_in = 1;
    cR_idx_in = eta_bar_idx_in + 1 : eta_bar_idx_in + magnetic_axis_data.top;
    sZ_idx_in = cR_idx_in(end) + 1 : cR_idx_in(end) + magnetic_axis_data.top;
    I_idx_in = sZ_idx_in(end) + 1 : sZ_idx_in(end) + numel(coilData.I);
    cc_idx_in = I_idx_in(end) + 1 : I_idx_in(end) + numel(coilData.coil_coeffs);
    idxs.eta_bar_idx = eta_bar_idx_in;
    idxs.cR_idx = cR_idx_in;
    idxs.sZ_idx = sZ_idx_in;
    idxs.I_idx = I_idx_in;
    idxs.cc_idx = cc_idx_in;
    

    dof = 1 + 2*magnetic_axis_data.top + numel(coilData.I)  + numel(coilData.coil_coeffs);
    x = zeros(1,dof);
    
    
    x(idxs.eta_bar_idx) = magnetic_axis_data.eta_bar ;
    x(idxs.cR_idx) = magnetic_axis_data.cR';
    x(idxs.sZ_idx) = magnetic_axis_data.sZ';
    x(idxs.I_idx)  = coilData.I'*gamma;
    x(idxs.cc_idx) = coilData.coil_coeffs(:)';
    
   
    
    inargs.filename_f = fopen(['./output/f_',datestr(now,'yyyymmdd_HH_MM_SS'),'.txt'],'w');
    inargs.idxs = idxs;

    
    
    tic
    inargs.maxiter = -1;
    options = optimoptions('fminunc','Display','iter-detailed','Algorithm','quasi-newton','SpecifyObjectiveGradient',true, 'TolFun', 0, 'TolX', 0, 'MaxIteration', Inf, 'OutputFcn',@(x,optimValues,state)output_function(x,optimValues,state, inargs));
    [x_out,fval,exitflag,output,grad] = fminunc(@(x) quasi(x, coilData, magnetic_axis_data, objective_data, idxs), x ,options);
    fprintf("\n\n done!\n");
    toc
 