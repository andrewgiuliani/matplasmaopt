function stop = output_function(x,optimValues,state, inargs)
        %% save previous iteration
        coil_idxs = [inargs.idxs.cc_idx, inargs.idxs.I_idx];
        axis_idxs = [inargs.idxs.sZ_idx, inargs.idxs.cR_idx, inargs.idxs.eta_bar_idx];
        fprintf(inargs.filename_f, '%.16e ',[optimValues.fval, optimValues.funccount, max(abs(optimValues.gradient(coil_idxs))), max(abs(optimValues.gradient(axis_idxs)))]);
        fprintf(inargs.filename_f,'\n');
        
        stop = false;

end

