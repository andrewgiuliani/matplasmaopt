classdef coilDataType < handle

    
    properties
        nfp
        ss
        coil_coeffs
        coil_points
        coil_field
        coil_tangents
        I
        C
        Nt
        PPP
        M
        fourierData
    end
    
    methods
        function obj = coilDataType(coil_coeffs_name, currents_name, in_nfp, in_ss, in_Nt)
            filename = coil_coeffs_name;
            filename_currents = currents_name;
            obj.coil_coeffs = dlmread(filename);
            obj.I = dlmread(filename_currents);
            
            obj.nfp = in_nfp;
            obj.ss = in_ss;
            
            % number of coils
            obj.C = size(obj.coil_coeffs, 2)/(3 * 2);
            % number of terms
            obj.Nt = size(obj.coil_coeffs, 1);
            if nargin > 4
                obj.Nt = in_Nt;
                obj.coil_coeffs = obj.coil_coeffs(1:obj.Nt, :);
            end
            
            obj.PPP = 20;% number of points per period
            obj.M = obj.PPP * obj.Nt;% number of discretization points for each coil

            
            
        
            %% precompute values of sines and cosines for the coils
            tt = linspace(0,2*pi,obj.M+1)';tt = tt(1:end-1);
            obj.fourierData.COS = zeros(obj.M,obj.Nt);
            obj.fourierData.SIN = zeros(obj.M,obj.Nt);
            obj.fourierData.dCOS = zeros(obj.M,obj.Nt);
            obj.fourierData.dSIN = zeros(obj.M,obj.Nt);

            for i = 1:obj.Nt   
                obj.fourierData.COS(:,i) = cos((i-1)*tt);
                obj.fourierData.SIN(:,i) = sin((i-1)*tt);
                obj.fourierData.dCOS(:,i) = -(i-1)*sin((i-1)*tt);
                obj.fourierData.dSIN(:,i) = (i-1)*cos((i-1)*tt);
            end
        end
        

        function [len] = get_clength(self)
            len = 0.;
            for k = 1: self.C
                len = len + sum(sqrt(self.coil_tangents{k}(:,1).^2+self.coil_tangents{k}(:,2).^2+self.coil_tangents{k}(:,3).^2));
            end
            len = len*(2*pi / self.M);
        end
        
        function [len] = get_clength_ind(self, k)
            len =  sum(sqrt(self.coil_tangents{k}(:,1).^2+self.coil_tangents{k}(:,2).^2+self.coil_tangents{k}(:,3).^2));
            len = len*(2*pi / self.M);
        end
        

        
        function [dlen_dc] = compute_dlen_dc(self)
            mag_T = zeros(self.M, self.C);
            for i = 1: self.C
                mag_T(:,i) = sqrt(self.coil_tangents{i}(:,1).^2+self.coil_tangents{i}(:,2).^2+self.coil_tangents{i}(:,3).^2);
            end
            
            
            dT_dc_dot_T = zeros(self.M,numel(self.coil_coeffs));
                
            for j = 1:numel(self.coil_coeffs)
                [k, kk] = ind2sub(size(self.coil_coeffs),j) ;

                i = ceil(kk/ 2/3);
                coord = mod(ceil(kk/2)-1,3)+1; %x,y,z

                if mod(kk,2) ~=0
                    dT_dc_dot_T(:,j) = self.fourierData.dSIN(:,k) .* self.coil_tangents{i}(:,coord)./mag_T(:,i);
                else
                    dT_dc_dot_T(:,j) = self.fourierData.dCOS(:,k) .* self.coil_tangents{i}(:,coord)./mag_T(:,i);
                end
            end
            
            dlen_dc = sum(dT_dc_dot_T,1) * 2 * pi / self.M;
        end

       function [cumul_len] = get_cumul_length_ind(self, k)
            dl =  (2*pi / self.M)*sqrt(self.coil_tangents{k}(:,1).^2+self.coil_tangents{k}(:,2).^2+self.coil_tangents{k}(:,3).^2);
            cumul_len = cumsum(dl);
       end
       
       
       
       
       
       
        function [z] = compute_position_ind(self, i, theta)              
            SIN = sin( (0:self.Nt-1) * theta );
            COS = cos( (0:self.Nt-1) * theta );
                
%             x = SIN*self.coil_coeffs(:,6*(i-1) + 1) + COS*self.coil_coeffs(:,6*(i-1) + 2);
%             y = SIN*self.coil_coeffs(:,6*(i-1) + 3) + COS*self.coil_coeffs(:,6*(i-1) + 4);
            z = SIN*self.coil_coeffs(:,6*(i-1) + 5) + COS*self.coil_coeffs(:,6*(i-1) + 6);
  
        end
        
        
        % z position must be 0 at theta = 0.
        function [pos_dcc] = compute_position_dcc(self)
            SIN =  sin( (0:self.Nt-1) * 0 );
            COS =  cos( (0:self.Nt-1) * 0 );
            pos_dcc = zeros(1,numel(self.coil_coeffs));
                
            for j = 1:numel(self.coil_coeffs)
                [k, kk] = ind2sub(size(self.coil_coeffs),j) ;
                i = ceil(kk/ 2/3);
                coord = mod(ceil(kk/2)-1,3)+1; %x,y,z
                
                if coord == 1 || coord == 2
                    continue
                end
                
                if mod(kk,2) ~=0
                    pos_dcc(1,j) = SIN(k) ;
                else
                    pos_dcc(1,j) = COS(k);
                end
            end
            
        end
       
       
        
        
        function cumul_length_term = get_cumul_length_term(self, c)
            tt = linspace(0,1.,self.M+1)';tt = tt(2:end);
            
           cumul_length_temp = self.get_cumul_length_ind(c);
           cumul_length = cumul_length_temp/cumul_length_temp(end);

           cumul_length_term = cumul_length - tt;

        end
        
        
        function cumul_length_dcc = compute_cumul_length_dcc(self)
            cumul_length_dcc = zeros(self.M, numel(self.coil_coeffs));
            
            mag_T = zeros(self.M, self.C);
            for i = 1: self.C
                mag_T(:,i) = sqrt(self.coil_tangents{i}(:,1).^2+self.coil_tangents{i}(:,2).^2+self.coil_tangents{i}(:,3).^2);
            end
            
            
            dT_dc = zeros(self.M,numel(self.coil_coeffs));
                
            for j = 1:numel(self.coil_coeffs)
                [k, kk] = ind2sub(size(self.coil_coeffs),j) ;

                i = ceil(kk/ 2/3);
                coord = mod(ceil(kk/2)-1,3)+1; %x,y,z

                if mod(kk,2) ~=0
                    dT_dc(:,j) = self.fourierData.dSIN(:,k) .* self.coil_tangents{i}(:,coord)./mag_T(:,i);
                else
                    dT_dc(:,j) = self.fourierData.dCOS(:,k) .* self.coil_tangents{i}(:,coord)./mag_T(:,i);
                end
            end
            
            L        = (2*pi / self.M) * sum(mag_T,1);
            dL_dc    = (2*pi / self.M) * sum(dT_dc,1);
            csl      = (2*pi / self.M) * cumsum(mag_T,1);
            d_csl_dc = (2*pi / self.M) * cumsum(dT_dc,1);
            
            for j = 1:numel(self.coil_coeffs)
                [k, kk] = ind2sub(size(self.coil_coeffs),j) ;

                i = ceil(kk/ 2/3);
                cumul_length_dcc(:,j) = (d_csl_dc(:,j) * L(i) - dL_dc(j)*csl(:,i) )/L(i)^2;
            end
            
        end
        
        
        
        
        
        function update(self)
            self.coil_points = zeros(self.M, 3*self.C);
            self.coil_field{1,self.C*self.nfp*(1+self.ss)} = [];
            self.coil_tangents{1,self.C*self.nfp*(1+self.ss)} = [];
            [self.coil_tangents{1,1:self.C*self.nfp*(1+self.ss)}] = deal(zeros(self.M,3));
            
            
            
            % update all coil points
            for i = 1:self.C
                self.coil_points(:,3*(i-1) + 1) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 1) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 2);
                self.coil_points(:,3*(i-1) + 2) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 3) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 4);
                self.coil_points(:,3*(i-1) + 3) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 5) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 6);
            end

            Ref = [1, 0 0.; 0, -1, 0; 0 0 -1]';
            % update all coil coil field values tangents the points
            for i = 1:self.C
                self.coil_field{i}(:,1) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 1) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 2);
                self.coil_field{i}(:,2) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 3) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 4);
                self.coil_field{i}(:,3) = self.fourierData.SIN*self.coil_coeffs(:,6*(i-1) + 5) + self.fourierData.COS*self.coil_coeffs(:,6*(i-1) + 6);
                
                self.coil_tangents{i}(:,1) = self.fourierData.dSIN*self.coil_coeffs(:,6*(i-1) + 1) + self.fourierData.dCOS*self.coil_coeffs(:,6*(i-1) + 2);
                self.coil_tangents{i}(:,2) = self.fourierData.dSIN*self.coil_coeffs(:,6*(i-1) + 3) + self.fourierData.dCOS*self.coil_coeffs(:,6*(i-1) + 4);
                self.coil_tangents{i}(:,3) = self.fourierData.dSIN*self.coil_coeffs(:,6*(i-1) + 5) + self.fourierData.dCOS*self.coil_coeffs(:,6*(i-1) + 6);

                if self.ss == 1
                    self.coil_field{self.nfp*self.C + i} = self.coil_field{i}*Ref;
                    self.coil_tangents{self.nfp*self.C + i} = self.coil_tangents{i}*Ref;
                end

                for t = 1:self.nfp-1
                    Rot = [cos(t * 2 * pi / self.nfp), - sin(t * 2 * pi / self.nfp) 0.; sin(t * 2 * pi / self.nfp), cos(t * 2 * pi / self.nfp), 0; 0 0 1]';
                    self.coil_field{t*self.C + i} = self.coil_field{i}*Rot;
                    self.coil_tangents{t*self.C + i} = self.coil_tangents{i}*Rot;

                    if self.ss == 1
                        self.coil_field{self.nfp*self.C +t*self.C + i} = self.coil_field{i}*Rot*Ref;
                        self.coil_tangents{self.nfp*self.C +t*self.C + i} = self.coil_tangents{i}*Rot*Ref;
                    end
                end 
            end

        end
        
        
        function plot_coils(self)
            self.update();
%             figure(1);
            grid on;view(3);hold on;
            colorstring = 'rbgmrc';
            for i = 1:self.C

                p = [self.coil_points(:,3*(i-1) + 1)';self.coil_points(:,3*(i-1) + 2)';self.coil_points(:,3*(i-1) + 3)'];
                plot3([p(1,:) p(1,1)], [p(2,:) p(2,1)], [p(3,:) p(3,1)], '-', 'Color', colorstring(i)); hold on;

                if self.ss == 1
                    R = [1, 0 0.; 0, -1, 0; 0 0 -1];
                    p = R * p;
                    plot3([p(1,:) p(1,1)], [p(2,:) p(2,1)], [p(3,:) p(3,1)], '-', 'Color', colorstring(i)); hold on;
                end
                for t = 1:self.nfp-1
                    R = [cos(t * 2 * pi / self.nfp), - sin(t * 2 * pi / self.nfp) 0.; sin(t * 2 * pi / self.nfp), cos(t * 2 * pi / self.nfp), 0; 0 0 1];
                    p = R * [self.coil_points(:,3*(i-1) + 1)';self.coil_points(:,3*(i-1) + 2)';self.coil_points(:,3*(i-1) + 3)'];
                    plot3([p(1,:) p(1,1)], [p(2,:) p(2,1)], [p(3,:) p(3,1)], '-', 'Color', colorstring(i)); hold on;
                    if self.ss == 1
                        R = [1, 0 0.; 0, -1, 0; 0 0 -1];
                        p = R * p;
                        plot3([p(1,:) p(1,1)], [p(2,:) p(2,1)], [p(3,:) p(3,1)], '-', 'Color', colorstring(i)); hold on;
                    end
                end

            end
            axis equal;
            
            xlabel('X', 'fontSize', 12);
            ylabel('Y', 'fontSize', 12);
            zlabel('Z', 'fontSize', 12);
        end
        
        function write_coils(self)
            outputFilename1 = ['./output/coilData_',datestr(now,'yyyymmdd_HH_MM_SS'),'.txt'];
            dlmwrite(outputFilename1, self.coil_coeffs, 'precision', 16);
            outputFilename1 = ['./output/I_',datestr(now,'yyyymmdd_HH_MM_SS'),'.txt'];
            dlmwrite(outputFilename1, self.I, 'precision', 16);
        end
        
    end
end

