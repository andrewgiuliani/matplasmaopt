classdef objectiveDataType < handle
    %OBJECTIVEDATATYPE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        current_g
        current_c_len
        current_c_len_ind
        current_ma_len
        current_data
        
        current_iota
        c_target
        ma_len_target
        iota_target
        current_grad
        current_gradnorm
    end
    
    methods
        function obj = objectiveDataType(in_c_target,in_ma_len_target,in_iota_target)
            obj.c_target = in_c_target;
            obj.ma_len_target = in_ma_len_target;
            obj.iota_target = in_iota_target;
        end
        
        function set_current_values(self,in_current_g,in_current_c, in_current_c_ind, in_current_ma_len,in_current_iota, in_current_gradnorm)
            self.current_g = in_current_g;
            self.current_c_len = in_current_c;
            self.current_c_len_ind = in_current_c_ind;
            self.current_ma_len = in_current_ma_len;
            self.current_iota = in_current_iota;
            self.current_gradnorm = in_current_gradnorm;
        end
    end
end

