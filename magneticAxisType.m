classdef magneticAxisType < handle

    properties   
        nfp
        
        % axis (R,Z) values
        RA
        ZA
        RAp
        ZAp
        RApp
        ZApp
        RAppp
        ZAppp
        
        dRA_dcos
        dRAp_dcos
        dRApp_dcos
        dRAppp_dcos
        
        dZA_dsin
        dZAp_dsin
        dZApp_dsin
        dZAppp_dsin
        
        sigma
        iota
        eta_bar
        curvature
        torsion
        B0
        sign_G
        sign_psi
        I2_over_B0
        abs_G_over_B0
        d_abs_G_over_B0_dRcos
        d_abs_G_over_B0_dZsin
        tau_kappa2
        
        s_theta
        s_thetatheta
        s_thetathetatheta
        
        X1c
        Y1s
        Y1c
        d_X1c
        d_Y1s
        d_Y1c
        
        d_Y1c_ds
        d2_Y1c_ds
        
        B_QS
        gradBx_QS
        gradBy_QS
        gradBz_QS
        
        B_coils
        gradBx_coils
        gradBy_coils
        gradBz_coils
        
        d_curvature_dtheta
        
        r_theta_cyl
        r_thetatheta_cyl
        r_thetathetatheta_cyl
        
        r_theta_cyl_tc
        r_thetatheta_cyl_tc
        
        kappa
        d_kappa_dtheta
        d_kappa_dtheta1
        d_kappa_dtheta2
        d_kappa_dtheta3
        
        % sine/cos values on axis
        cR
        sZ
        D
        Dtilde
        zeta
        nzeta
        S
        C
        SIN
        COS
        dSIN
        dCOS
        top 

        
        % on axis magnetic field values and sensitivities
        B
        B_mag
        BR
        Bphi
        BZ
        t_vec 
        n_vec
        b_vec
        dX_dPhi % Tx
        dY_dPhi % Ty
        dZ_dPhi % TZ
        dTx_dr
        dTy_dr
        dTz_dr
        dTx_dz
        dTy_dz
        dTz_dz
        Ta_mag
        dB_R_dc
        dB_phi_dc
        dB_Z_dc
        dR_dr
        dRp_dr
        
        d_Ta_mag_dRcos
        d_Ta_mag_dZsin
        
        d_s_theta_dZsin
        d_s_theta_dRcos
        d_s_thetatheta_dRcos
        d_s_thetatheta_dZsin
        d_s_thetathetatheta_dRcos
        d_s_thetathetatheta_dZsin
        
        d_r_theta1_cyl_dRcos
        d_r_theta2_cyl_dRcos
        d_r_theta3_cyl_dRcos

        d_r_theta1_cyl_dZsin
        d_r_theta2_cyl_dZsin
        d_r_theta3_cyl_dZsin
        
        d_r_theta1_cyl_tc_dRcos
        d_r_theta2_cyl_tc_dRcos
        d_r_theta3_cyl_tc_dRcos

        d_r_theta1_cyl_tc_dZsin
        d_r_theta2_cyl_tc_dZsin
        d_r_theta3_cyl_tc_dZsin

        d_r_thetatheta1_cyl_dRcos
        d_r_thetatheta2_cyl_dRcos
        d_r_thetatheta3_cyl_dRcos

        d_r_thetatheta1_cyl_dZsin
        d_r_thetatheta2_cyl_dZsin
        d_r_thetatheta3_cyl_dZsin

        d_r_thetatheta1_cyl_tc_dRcos
        d_r_thetatheta2_cyl_tc_dRcos
        d_r_thetatheta3_cyl_tc_dRcos

        d_r_thetatheta1_cyl_tc_dZsin
        d_r_thetatheta2_cyl_tc_dZsin
        d_r_thetatheta3_cyl_tc_dZsin

        d_r_thetathetatheta1_cyl_dRcos
        d_r_thetathetatheta2_cyl_dRcos
        d_r_thetathetatheta3_cyl_dRcos

        d_r_thetathetatheta1_cyl_dZsin
        d_r_thetathetatheta2_cyl_dZsin
        d_r_thetathetatheta3_cyl_dZsin
        
        
        
        d_curvature_dRcos
        d_curvature_dZsin
        
        d_kappaR_dRcos
        d_kappaPhi_dRcos
        d_kappaZ_dRcos
        d_kappaR_dZsin
        d_kappaPhi_dZsin
        d_kappaZ_dZsin
        
        
        dt1_dRcos
        dt2_dRcos
        dt3_dRcos
        dt1_dZsin
        dt2_dZsin
        dt3_dZsin
        
        dn1_dRcos
        dn2_dRcos
        dn3_dRcos
        dn1_dZsin
        dn2_dZsin
        dn3_dZsin
        
        db1_dRcos
        db2_dRcos
        db3_dRcos
        db1_dZsin
        db2_dZsin
        db3_dZsin
        
        d_kappa_dtheta1_dRcos
        d_kappa_dtheta2_dRcos
        d_kappa_dtheta3_dRcos

        d_kappa_dtheta1_dZsin
        d_kappa_dtheta2_dZsin
        d_kappa_dtheta3_dZsin

        d_curvature_dtheta_dRcos
        d_curvature_dtheta_dZsin

        d_torsion_dRcos
        d_torsion_dZsin
        
        torsion_numerator
        torsion_denominator
        
        d_tau_kappa2_dRcos
        d_tau_kappa2_dZsin                             
        d_abs_G0B0taukappa2_dRcos
        d_abs_G0B0taukappa2_dZsin
        
        
        
        dBx_dRcos
        dBy_dRcos
        dBz_dRcos
        dBx_dZsin
        dBy_dZsin
        dBz_dZsin
        d_X1c_dRcos
        d_X1c_dZsin
        d_Y1s_dRcos
        d_Y1s_dZsin
        d_Y1c_dRcos
        d_Y1c_dZsin
        d2_X1c_dRcos
        d2_X1c_dZsin
        d2_Y1s_dRcos
        d2_Y1s_dZsin
        d2_Y1c_dRcos
        d2_Y1c_dZsin
        
        d_B1x_dRcos
        d_B1y_dRcos
        d_B1z_dRcos
        d_B2x_dRcos
        d_B2y_dRcos
        d_B2z_dRcos
        d_B3x_dRcos
        d_B3y_dRcos
        d_B3z_dRcos

        d_B1x_dZsin
        d_B1y_dZsin
        d_B1z_dZsin
        d_B2x_dZsin
        d_B2y_dZsin
        d_B2z_dZsin
        d_B3x_dZsin
        d_B3y_dZsin
        d_B3z_dZsin
        
        dBx_coils_dc
        dBy_coils_dc
        dBz_coils_dc
        
        dBx_coils_dRcos
        dBy_coils_dRcos
        dBz_coils_dRcos
        
        
        
        dBx_coils_dZsin
        dBy_coils_dZsin
        dBz_coils_dZsin
        
        dBx_coils_dI
        dBy_coils_dI
        dBz_coils_dI
        
        dB1x_dI
        dB1y_dI
        dB1z_dI
        dB2x_dI
        dB2y_dI
        dB2z_dI
        dB3x_dI
        dB3y_dI
        dB3z_dI
        
        dB1x_dRcos
        dB1y_dRcos
        dB1z_dRcos
        dB2x_dRcos
        dB2y_dRcos
        dB2z_dRcos
        dB3x_dRcos
        dB3y_dRcos
        dB3z_dRcos



        dB1x_dZsin
        dB1y_dZsin
        dB1z_dZsin
        dB2x_dZsin
        dB2y_dZsin
        dB2z_dZsin
        dB3x_dZsin
        dB3y_dZsin
        dB3z_dZsin
        
        B1x_coils_dc
        B1y_coils_dc
        B1z_coils_dc
        
        B2x_coils_dc
        B2y_coils_dc
        B2z_coils_dc
        
        B3x_coils_dc
        B3y_coils_dc
        B3z_coils_dc

        
        d_B1x_ds
        d_B1y_ds
        d_B1z_ds

        d_B2x_ds
        d_B2y_ds
        d_B2z_ds

        d_B3x_ds
        d_B3y_ds
        d_B3z_ds
        
        d_gradB1_diota
        d_gradB2_diota
        d_gradB3_diota
        
        d_gradB1_detabar
        d_gradB2_detabar
        d_gradB3_detabar
        
        d_X1c_detabar
        d_Y1s_detabar 
        d_Y1c_detabar
        
        d2_X1c_detabar
        d2_Y1s_detabar 
        d2_Y1c_detabar
        
        svp
        sg 
    end
    
    methods(Static)
        function [R,Z] = get_coords(self,phi)
            num_modes = size(self.cR,1);
            
            R = 0.;
            Z = 0.;
            for j = 1:num_modes
                R = R + self.cR(j)*cos(self.nfp * phi * (j-1));
                Z = Z + self.sZ(j)*sin(self.nfp * phi * (j-1));
            end
        end
        function [dR_dphi,dZ_dphi] = get_dcoords(self,phi)
            num_modes = size(self.cR,1);
            
            dR_dphi = 0.;
            dZ_dphi = 0.;
            for j = 1:num_modes
                dR_dphi = dR_dphi -self.cR(j)*(j-1)*self.nfp * sin(self.nfp * phi * (j-1));
                dZ_dphi = dZ_dphi +self.sZ(j)*(j-1)*self.nfp * cos(self.nfp * phi * (j-1));
            end
        end
    end
    methods  
        function obj = magneticAxisType(coilData, in_cR, in_sZ, in_top)
            obj.nfp = coilData.nfp;
            %% differentiation matrices
            obj.nzeta  = coilData.M;
            if mod(obj.nzeta,2)==0
                obj.nzeta = obj.nzeta+1;
            end
            scheme = 20;
            [obj.zeta, ~, obj.D,  ~] = sfincs_uniformDiffMatrices(obj.nzeta,  0, 2*pi/obj.nfp, scheme);
            
            obj.cR = dlmread(in_cR);
            obj.sZ = dlmread(in_sZ);
            
            obj.cR = obj.cR(1:in_top);
            obj.sZ = obj.sZ(1:in_top);
            
            obj.top = in_top;
            

            
            %% precompute values of sines and cosines on the magnetic axis
            tt = linspace(0,2*pi/obj.nfp,obj.nzeta+1)';tt = tt(1:end-1);
            obj.COS = zeros(obj.nzeta,obj.top);
            obj.SIN = zeros(obj.nzeta,obj.top);
            obj.dCOS = zeros(obj.nzeta,obj.top);
            obj.dSIN = zeros(obj.nzeta,obj.top);

            for i = 1:obj.top
                obj.COS(:,i) = cos(obj.nfp*(i-1)*tt);
                obj.SIN(:,i) = sin(obj.nfp*(i-1)*tt);
                obj.dCOS(:,i) = -obj.nfp*(i-1)*sin(obj.nfp*(i-1)*tt);
                obj.dSIN(:,i) = obj.nfp*(i-1)*cos(obj.nfp*(i-1)*tt);
            end
            
            obj.B = zeros(obj.nzeta,3);
            obj.B_mag = zeros(obj.nzeta,1);
            obj.BR = zeros(obj.nzeta,3);
            obj.Bphi = zeros(obj.nzeta,3);
            obj.BZ = zeros(obj.nzeta,3);
            obj.dX_dPhi = zeros(obj.nzeta,1);% Tx
            obj.dY_dPhi = zeros(obj.nzeta,1);% Ty
            obj.dZ_dPhi = zeros(obj.nzeta,1);% TZ
            obj.Ta_mag = zeros(obj.nzeta,1);
            
            obj.dTx_dr = zeros(obj.nzeta,obj.nzeta);
            obj.dTy_dr = zeros(obj.nzeta,obj.nzeta);
            obj.dTz_dr = zeros(obj.nzeta,obj.nzeta);
            obj.dTx_dz = zeros(obj.nzeta,obj.nzeta);
            obj.dTy_dz = zeros(obj.nzeta,obj.nzeta);
            obj.dTz_dz = zeros(obj.nzeta,obj.nzeta);

            obj.dB_R_dc = zeros(obj.nzeta, numel(coilData.coil_coeffs));
            obj.dB_phi_dc = zeros(obj.nzeta, numel(coilData.coil_coeffs));
            obj.dB_Z_dc = zeros(obj.nzeta, numel(coilData.coil_coeffs));
            
            obj.dR_dr = zeros(obj.nzeta, obj.nzeta);
            obj.dRp_dr = zeros(obj.nzeta, obj.nzeta);
            
            obj.sigma = zeros(obj.nzeta,1);
            obj.iota = 0.;
            obj.eta_bar = 1;
            obj.B0 = 1.;
            
           
            obj.RA = zeros(obj.nzeta,1);
            obj.ZA = zeros(obj.nzeta,1);
            obj.RAp = zeros(obj.nzeta,1);
            obj.ZAp = zeros(obj.nzeta,1);
            obj.RApp = zeros(obj.nzeta,1);
            obj.ZApp = zeros(obj.nzeta,1);
            obj.RAppp = zeros(obj.nzeta,1);
            obj.ZAppp = zeros(obj.nzeta,1);

            obj.dRA_dcos = zeros(obj.nzeta,obj.top);
            obj.dRAp_dcos = zeros(obj.nzeta,obj.top);
            obj.dRApp_dcos = zeros(obj.nzeta,obj.top);
            obj.dRAppp_dcos = zeros(obj.nzeta,obj.top);


            obj.dZA_dsin = zeros(obj.nzeta,obj.top);
            obj.dZAp_dsin = zeros(obj.nzeta,obj.top);
            obj.dZApp_dsin = zeros(obj.nzeta,obj.top);
            obj.dZAppp_dsin = zeros(obj.nzeta,obj.top);
            for imn = 1:obj.top
                n = obj.nfp*(imn-1);
                angle = n * obj.zeta;
                sinangle = sin(angle);
                cosangle = cos(angle);
                obj.dRA_dcos(:,imn) = cosangle;
                obj.dZA_dsin(:,imn) = sinangle;

                obj.dRAp_dcos(:,imn) = -n*sinangle;
                obj.dZAp_dsin(:,imn) =  n*cosangle;

                obj.dRApp_dcos(:,imn) = - n*n*cosangle;
                obj.dZApp_dsin(:,imn) = - n*n*sinangle;

                obj.dRAppp_dcos(:,imn) =   n*n*n*sinangle;
                obj.dZAppp_dsin(:,imn) = - n*n*n*cosangle;
            end
            
            
            obj.dBx_coils_dc =  zeros(obj.nzeta, numel(coilData.coil_coeffs));
            obj.dBy_coils_dc =  zeros(obj.nzeta, numel(coilData.coil_coeffs));
            obj.dBz_coils_dc =  zeros(obj.nzeta, numel(coilData.coil_coeffs));
            
            
            obj.B1x_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B1y_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B1z_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));

            obj.B2x_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B2y_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B2z_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));

            obj.B3x_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B3y_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));
            obj.B3z_coils_dc = zeros(obj.nzeta,  numel(coilData.coil_coeffs));

            obj.dBx_coils_dRcos = zeros(obj.nzeta, obj.top);
            obj.dBy_coils_dRcos = zeros(obj.nzeta, obj.top);
            obj.dBz_coils_dRcos = zeros(obj.nzeta, obj.top);



            obj.dBx_coils_dZsin = zeros(obj.nzeta, obj.top);
            obj.dBy_coils_dZsin = zeros(obj.nzeta, obj.top);
            obj.dBz_coils_dZsin = zeros(obj.nzeta, obj.top);
            
            obj.dBx_coils_dI = zeros(obj.nzeta,size(coilData.I,1));
            obj.dBy_coils_dI = zeros(obj.nzeta,size(coilData.I,1));
            obj.dBz_coils_dI = zeros(obj.nzeta,size(coilData.I,1));
            
            obj.dB1x_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB1y_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB1z_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB2x_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB2y_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB2z_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB3x_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB3y_dI = zeros(obj.nzeta, size(coilData.I,1));
            obj.dB3z_dI = zeros(obj.nzeta, size(coilData.I,1));
            
            obj.dB1x_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB1y_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB1z_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB2x_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB2y_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB2z_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB3x_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB3y_dRcos = zeros(obj.nzeta, obj.top);
            obj.dB3z_dRcos = zeros(obj.nzeta, obj.top);



            obj.dB1x_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB1y_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB1z_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB2x_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB2y_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB2z_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB3x_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB3y_dZsin = zeros(obj.nzeta, obj.top);
            obj.dB3z_dZsin = zeros(obj.nzeta, obj.top);
            
            obj.d_gradB1_detabar = zeros(obj.nzeta,3);
            obj.d_gradB2_detabar = zeros(obj.nzeta,3);
            obj.d_gradB3_detabar = zeros(obj.nzeta,3);

            
            obj.svp = 1;
            obj.sg  = 1;
        end
        
       


        
        function update(self)    
            
            d2X_dPhi_dRcos = self.dRAp_dcos .* cos(self.zeta) - self.dRA_dcos .* sin(self.zeta);
            d2Y_dPhi_dRcos = self.dRAp_dcos .* sin(self.zeta) + self.dRA_dcos .* cos(self.zeta);
            d2Z_dPhi_dRcos = 0.;
            self.d_Ta_mag_dRcos = (self.dX_dPhi .* d2X_dPhi_dRcos + self.dY_dPhi .* d2Y_dPhi_dRcos + self.dZ_dPhi .* d2Z_dPhi_dRcos)./self.Ta_mag;
            
            d2X_dPhi_dZsin = 0.;
            d2Y_dPhi_dZsin = 0.;
            d2Z_dPhi_dZsin = self.dZAp_dsin;
            self.d_Ta_mag_dZsin = (self.dX_dPhi .* d2X_dPhi_dZsin + self.dY_dPhi .* d2Y_dPhi_dZsin + self.dZ_dPhi .* d2Z_dPhi_dZsin)./self.Ta_mag;
            
            
            
            % compute all partial derivatives wrt to ca on the magnetic axis
            % assumes that compute_f has already been called
            self.d_s_theta_dZsin = zeros(self.nzeta,numel(self.sZ));
            self.d_s_theta_dRcos = zeros(self.nzeta,numel(self.cR));

            self.d_s_thetatheta_dRcos = zeros(self.nzeta,numel(self.cR));
            self.d_s_thetatheta_dZsin = zeros(self.nzeta,numel(self.sZ));


            for imn = 1:numel(self.cR)
                self.d_s_theta_dRcos(:,imn) = (self.RA.*self.dRA_dcos(:,imn) + self.RAp.*self.dRAp_dcos(:,imn) )./ self.s_theta;
                self.d_s_theta_dZsin(:,imn) = self.ZAp.*self.dZAp_dsin(:,imn)./ self.s_theta;
            end

            for imn = 1:numel(self.cR)
                self.d_s_thetatheta_dRcos(:,imn) = self.RA.*(self.dRAp_dcos(:,imn)./self.s_theta - self.RAp.*self.d_s_theta_dRcos(:,imn)./self.s_theta.^2)...
                                             +self.RAp.*((self.dRA_dcos(:,imn)+self.dRApp_dcos(:,imn))./self.s_theta - self.RApp.*self.d_s_theta_dRcos(:,imn)./self.s_theta.^2)...
                                             +self.ZAp.*( - self.ZApp.*self.d_s_theta_dRcos(:,imn)./self.s_theta.^2 ) ...
                                             +(self.dRAp_dcos(:,imn).*self.RApp )./self.s_theta;                               
                self.d_s_thetatheta_dZsin(:,imn) = self.RA.*( - self.RAp.*self.d_s_theta_dZsin(:,imn)./self.s_theta.^2)...
                                             +self.RAp.*( - self.RApp.*self.d_s_theta_dZsin(:,imn)./self.s_theta.^2)...
                                             +self.ZAp.*( self.dZApp_dsin(:,imn)./self.s_theta - self.ZApp.*self.d_s_theta_dZsin(:,imn)./self.s_theta.^2 )...
                                             +( self.dZAp_dsin(:,imn).*self.ZApp )./self.s_theta;
            end
                                                
            self.d_s_thetathetatheta_dRcos =  compute_d_stt_dRc(self.RA, self.RAp, self.RApp, self.RAppp, self.ZA, self.ZAp, self.ZApp, self.ZAppp,...
                                                                self.dRA_dcos,self.dRAp_dcos,self.dRApp_dcos,self.dRAppp_dcos,...
                                                                self.s_theta, self.s_thetatheta, self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,...
                                                                self.zeta, self.cR);                                                                                    
            self.d_s_thetathetatheta_dZsin =  compute_d_stt_dZc(self.RA, self.RAp, self.RApp, self.RAppp, self.ZA, self.ZAp, self.ZApp, self.ZAppp,...
                                                                self.dZA_dsin, self.dZAp_dsin, self.dZApp_dsin, self.dZAppp_dsin,...
                                                                self.s_theta, self.s_thetatheta, self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,...
                                                                self.zeta, self.sZ);   


            self.d_abs_G_over_B0_dRcos = mean(self.d_s_theta_dRcos,1)/self.B0;
            self.d_abs_G_over_B0_dZsin = mean(self.d_s_theta_dZsin,1)/self.B0;
   
            self.d_r_theta1_cyl_dRcos = self.dRAp_dcos;
            self.d_r_theta2_cyl_dRcos = self.dRA_dcos;
            self.d_r_theta3_cyl_dRcos = zeros(self.nzeta,numel(self.cR));

            self.d_r_theta1_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_theta2_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_theta3_cyl_dZsin = self.dZAp_dsin;

            self.d_r_theta1_cyl_tc_dRcos = self.dRApp_dcos;
            self.d_r_theta2_cyl_tc_dRcos = self.dRAp_dcos;
            self.d_r_theta3_cyl_tc_dRcos = zeros(self.nzeta,numel(self.cR));


            self.d_r_theta1_cyl_tc_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_theta2_cyl_tc_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_theta3_cyl_tc_dZsin = self.dZApp_dsin;

            self.d_r_thetatheta1_cyl_dRcos = self.dRApp_dcos - self.dRA_dcos;
            self.d_r_thetatheta2_cyl_dRcos = 2*self.dRAp_dcos;
            self.d_r_thetatheta3_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));


            self.d_r_thetatheta1_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetatheta2_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetatheta3_cyl_dZsin = self.dZApp_dsin;


            self.d_r_thetatheta1_cyl_tc_dRcos = self.dRAppp_dcos - self.dRAp_dcos;
            self.d_r_thetatheta2_cyl_tc_dRcos = 2*self.dRApp_dcos;
            self.d_r_thetatheta3_cyl_tc_dRcos =  zeros(self.nzeta,numel(self.cR));

            self.d_r_thetatheta1_cyl_tc_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetatheta2_cyl_tc_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetatheta3_cyl_tc_dZsin = self.dZAppp_dsin;


            self.d_r_thetathetatheta1_cyl_dRcos = self.dRAppp_dcos-3*self.dRAp_dcos;
            self.d_r_thetathetatheta2_cyl_dRcos = 3.*self.dRApp_dcos - self.dRA_dcos;
            self.d_r_thetathetatheta3_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));

            self.d_r_thetathetatheta1_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetathetatheta2_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.d_r_thetathetatheta3_cyl_dZsin = self.dZAppp_dsin;

            
            
            
            self.d_kappaR_dRcos = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                            self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,...
                                            self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1), ...
                                            self.d_r_theta1_cyl_dRcos, self.d_r_thetatheta1_cyl_dRcos,...
                                            self.zeta, self.cR);
            self.d_kappaPhi_dRcos = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                              self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,...
                                              self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2), ...
                                              self.d_r_theta2_cyl_dRcos, self.d_r_thetatheta2_cyl_dRcos,...
                                            self.zeta, self.cR);
            self.d_kappaZ_dRcos = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                            self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,...
                                            self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3), ...
                                            self.d_r_theta3_cyl_dRcos, self.d_r_thetatheta3_cyl_dRcos,...
                                            self.zeta, self.cR);                         
 



                          
            self.d_kappaR_dZsin = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                            self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,...
                                            self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1), ...
                                            self.d_r_theta1_cyl_dZsin, self.d_r_thetatheta1_cyl_dZsin,...
                                            self.zeta, self.sZ);
            self.d_kappaPhi_dZsin = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                              self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,...
                                              self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2), ...
                                              self.d_r_theta2_cyl_dZsin, self.d_r_thetatheta2_cyl_dZsin,...
                                            self.zeta, self.sZ);
            self.d_kappaZ_dZsin = compute_d_k_dc(self.s_theta, self.s_thetatheta,  ...
                                            self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,...
                                            self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3), ...
                                            self.d_r_theta3_cyl_dZsin, self.d_r_thetatheta3_cyl_dZsin,...
                                            self.zeta, self.sZ);                                       







            self.d_curvature_dRcos = (self.kappa(:,1) .* self.d_kappaR_dRcos + self.kappa(:,2) .* self.d_kappaPhi_dRcos + self.kappa(:,3) .* self.d_kappaZ_dRcos)./self.curvature;
            self.d_curvature_dZsin = (self.kappa(:,1) .* self.d_kappaR_dZsin + self.kappa(:,2) .* self.d_kappaPhi_dZsin + self.kappa(:,3) .* self.d_kappaZ_dZsin)./self.curvature; 


            dt1_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));
            dt2_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));
            dt3_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));


            dt1_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            dt2_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            dt3_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));


            self.dt1_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.dt2_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.dt3_dRcos =  zeros(self.nzeta,numel(self.cR));


            self.dt1_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.dt2_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.dt3_dZsin =  zeros(self.nzeta,numel(self.sZ));



            for imn = 1:numel(self.cR)
                dt1_cyl_dRcos(:,imn) = (self.d_r_theta1_cyl_dRcos(:,imn).*self.s_theta - self.d_s_theta_dRcos(:,imn).*self.RAp)./self.s_theta.^2;
                dt2_cyl_dRcos(:,imn) = (self.d_r_theta2_cyl_dRcos(:,imn).*self.s_theta - self.d_s_theta_dRcos(:,imn).*self.RA) ./self.s_theta.^2;
                dt3_cyl_dRcos(:,imn) = (self.d_r_theta3_cyl_dRcos(:,imn).*self.s_theta - self.d_s_theta_dRcos(:,imn).*self.ZAp)./self.s_theta.^2;

                dt1_cyl_dZsin(:,imn) = (self.d_r_theta1_cyl_dZsin(:,imn).*self.s_theta - self.d_s_theta_dZsin(:,imn).*self.RAp)./self.s_theta.^2;
                dt2_cyl_dZsin(:,imn) = (self.d_r_theta2_cyl_dZsin(:,imn).*self.s_theta - self.d_s_theta_dZsin(:,imn).*self.RA) ./self.s_theta.^2;
                dt3_cyl_dZsin(:,imn) = (self.d_r_theta3_cyl_dZsin(:,imn).*self.s_theta - self.d_s_theta_dZsin(:,imn).*self.ZAp)./self.s_theta.^2;
            end

            for imn = 1:numel(self.cR)
                self.dt1_dRcos(:,imn) = cos(self.zeta).*dt1_cyl_dRcos(:,imn)-sin(self.zeta).*dt2_cyl_dRcos(:,imn);
                self.dt2_dRcos(:,imn) = sin(self.zeta).*dt1_cyl_dRcos(:,imn)+cos(self.zeta).*dt2_cyl_dRcos(:,imn);
                self.dt3_dRcos(:,imn) = dt3_cyl_dRcos(:,imn);


                self.dt1_dZsin(:,imn) = cos(self.zeta).*dt1_cyl_dZsin(:,imn)-sin(self.zeta).*dt2_cyl_dZsin(:,imn);
                self.dt2_dZsin(:,imn) = sin(self.zeta).*dt1_cyl_dZsin(:,imn)+cos(self.zeta).*dt2_cyl_dZsin(:,imn);
                self.dt3_dZsin(:,imn) = dt3_cyl_dZsin(:,imn);
            end
        
            dn1_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));
            dn2_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));
            dn3_cyl_dRcos =  zeros(self.nzeta,numel(self.cR));



            dn1_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            dn2_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));
            dn3_cyl_dZsin =  zeros(self.nzeta,numel(self.sZ));


            self.dn1_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.dn2_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.dn3_dRcos =  zeros(self.nzeta,numel(self.cR));


            self.dn1_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.dn2_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.dn3_dZsin =  zeros(self.nzeta,numel(self.sZ));


            for imn = 1:numel(self.cR)
                dn1_cyl_dRcos(:,imn) = (self.d_kappaR_dRcos(:,imn).*self.curvature - self.d_curvature_dRcos(:,imn).*self.kappa(:,1))./self.curvature.^2;
                dn2_cyl_dRcos(:,imn)=(self.d_kappaPhi_dRcos(:,imn).*self.curvature - self.d_curvature_dRcos(:,imn).*self.kappa(:,2))./self.curvature.^2;
                dn3_cyl_dRcos(:,imn) = (self.d_kappaZ_dRcos(:,imn).*self.curvature - self.d_curvature_dRcos(:,imn).*self.kappa(:,3))./self.curvature.^2;

                dn1_cyl_dZsin(:,imn) = (self.d_kappaR_dZsin(:,imn).*self.curvature - self.d_curvature_dZsin(:,imn).*self.kappa(:,1))./self.curvature.^2;
                dn2_cyl_dZsin(:,imn)=(self.d_kappaPhi_dZsin(:,imn).*self.curvature - self.d_curvature_dZsin(:,imn).*self.kappa(:,2)) ./self.curvature.^2;
                dn3_cyl_dZsin(:,imn) = (self.d_kappaZ_dZsin(:,imn).*self.curvature - self.d_curvature_dZsin(:,imn).*self.kappa(:,3))./self.curvature.^2;
            end

            for imn = 1:numel(self.cR)
                self.dn1_dRcos(:,imn) = cos(self.zeta).*dn1_cyl_dRcos(:,imn)-sin(self.zeta).*dn2_cyl_dRcos(:,imn);
                self.dn2_dRcos(:,imn) = sin(self.zeta).*dn1_cyl_dRcos(:,imn)+cos(self.zeta).*dn2_cyl_dRcos(:,imn);
                self.dn3_dRcos(:,imn) = dn3_cyl_dRcos(:,imn);


                self.dn1_dZsin(:,imn) = cos(self.zeta).*dn1_cyl_dZsin(:,imn)-sin(self.zeta).*dn2_cyl_dZsin(:,imn);
                self.dn2_dZsin(:,imn) = sin(self.zeta).*dn1_cyl_dZsin(:,imn)+cos(self.zeta).*dn2_cyl_dZsin(:,imn);
                self.dn3_dZsin(:,imn) = dn3_cyl_dZsin(:,imn);
            end

            self.db1_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.db2_dRcos =  zeros(self.nzeta,numel(self.cR));
            self.db3_dRcos =  zeros(self.nzeta,numel(self.cR));


            self.db1_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.db2_dZsin =  zeros(self.nzeta,numel(self.sZ));
            self.db3_dZsin =  zeros(self.nzeta,numel(self.sZ));


            for imn = 1:numel(self.cR)
                    out =  cross([self.dt1_dRcos(:,imn), self.dt2_dRcos(:,imn), self.dt3_dRcos(:,imn)],self.n_vec)+cross(self.t_vec,[self.dn1_dRcos(:,imn), self.dn2_dRcos(:,imn), self.dn3_dRcos(:,imn)]);
                    self.db1_dRcos(:,imn) = out(:,1);
                    self.db2_dRcos(:,imn) = out(:,2);
                    self.db3_dRcos(:,imn) = out(:,3);


                    out =  cross([self.dt1_dZsin(:,imn), self.dt2_dZsin(:,imn), self.dt3_dZsin(:,imn)],self.n_vec)+cross(self.t_vec,[self.dn1_dZsin(:,imn), self.dn2_dZsin(:,imn), self.dn3_dZsin(:,imn)]);
                    self.db1_dZsin(:,imn) = out(:,1);
                    self.db2_dZsin(:,imn) = out(:,2);
                    self.db3_dZsin(:,imn) = out(:,3);
            end
            
            
        self.d_kappa_dtheta1_dRcos = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,self.d_s_thetathetatheta_dRcos,...
                                                  self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1), self.r_theta_cyl_tc(:,1), self.r_thetatheta_cyl_tc(:,1),...
                                                  self.d_r_theta1_cyl_dRcos, self.d_r_thetatheta1_cyl_dRcos, self.d_r_theta1_cyl_tc_dRcos,self.d_r_thetatheta1_cyl_tc_dRcos,...
                                                  self.zeta, self.cR);
        self.d_kappa_dtheta2_dRcos = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,self.d_s_thetathetatheta_dRcos,...
                                                  self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2), self.r_theta_cyl_tc(:,2), self.r_thetatheta_cyl_tc(:,2),...
                                                  self.d_r_theta2_cyl_dRcos, self.d_r_thetatheta2_cyl_dRcos, self.d_r_theta2_cyl_tc_dRcos,self.d_r_thetatheta2_cyl_tc_dRcos,...
                                                  self.zeta, self.cR);
        self.d_kappa_dtheta3_dRcos = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dRcos, self.d_s_thetatheta_dRcos,self.d_s_thetathetatheta_dRcos,...
                                                  self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3), self.r_theta_cyl_tc(:,3), self.r_thetatheta_cyl_tc(:,3),...
                                                  self.d_r_theta3_cyl_dRcos, self.d_r_thetatheta3_cyl_dRcos, self.d_r_theta3_cyl_tc_dRcos,self.d_r_thetatheta3_cyl_tc_dRcos,...
                                                  self.zeta, self.cR);

                                              
                                              
                                              
        self.d_kappa_dtheta1_dZsin = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,self.d_s_thetathetatheta_dZsin,...
                                                  self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1), self.r_theta_cyl_tc(:,1), self.r_thetatheta_cyl_tc(:,1),...
                                                  self.d_r_theta1_cyl_dZsin, self.d_r_thetatheta1_cyl_dZsin, self.d_r_theta1_cyl_tc_dZsin,self.d_r_thetatheta1_cyl_tc_dZsin,...
                                                  self.zeta, self.sZ);
        self.d_kappa_dtheta2_dZsin = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,self.d_s_thetathetatheta_dZsin,...
                                                  self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2), self.r_theta_cyl_tc(:,2), self.r_thetatheta_cyl_tc(:,2),...
                                                  self.d_r_theta2_cyl_dZsin, self.d_r_thetatheta2_cyl_dZsin, self.d_r_theta2_cyl_tc_dZsin,self.d_r_thetatheta2_cyl_tc_dZsin,...
                                                  self.zeta, self.sZ);
        self.d_kappa_dtheta3_dZsin = d_kappa_dtheta_dc(self.s_theta, self.s_thetatheta, self.s_thetathetatheta, ...
                                                  self.d_s_theta_dZsin, self.d_s_thetatheta_dZsin,self.d_s_thetathetatheta_dZsin,...
                                                  self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3), self.r_theta_cyl_tc(:,3), self.r_thetatheta_cyl_tc(:,3),...
                                                  self.d_r_theta3_cyl_dZsin, self.d_r_thetatheta3_cyl_dZsin, self.d_r_theta3_cyl_tc_dZsin,self.d_r_thetatheta3_cyl_tc_dZsin,...
                                                  self.zeta, self.sZ);
            
        self.d_curvature_dtheta_dRcos = compute_dk_dt_dc(self.kappa(:,1), self.kappa(:,2), self.kappa(:,3),...
                                                        self.d_kappa_dtheta(:,1), self.d_kappa_dtheta(:,2), self.d_kappa_dtheta(:,3),...
                                                        self.d_kappaR_dRcos, self.d_kappaPhi_dRcos, self.d_kappaZ_dRcos,...
                                                        self.d_kappa_dtheta1_dRcos, self.d_kappa_dtheta2_dRcos, self.d_kappa_dtheta3_dRcos,...
                                                        self.curvature,self.d_curvature_dRcos,...
                                                        self.zeta, self.cR);


        self.d_curvature_dtheta_dZsin = compute_dk_dt_dc(self.kappa(:,1), self.kappa(:,2), self.kappa(:,3),...
                                                        self.d_kappa_dtheta(:,1), self.d_kappa_dtheta(:,2), self.d_kappa_dtheta(:,3),...
                                                        self.d_kappaR_dZsin, self.d_kappaPhi_dZsin, self.d_kappaZ_dZsin,...
                                                        self.d_kappa_dtheta1_dZsin, self.d_kappa_dtheta2_dZsin, self.d_kappa_dtheta3_dZsin,...
                                                        self.curvature,self.d_curvature_dZsin,...
                                                        self.zeta, self.sZ);        
        
        
        


            self.d_torsion_dRcos = compute_torsion_dc(self.torsion_numerator,self.torsion_denominator, ...
                                                     self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1),self.r_thetathetatheta_cyl(:,1),self.d_r_theta1_cyl_dRcos, self.d_r_thetatheta1_cyl_dRcos, self.d_r_thetathetatheta1_cyl_dRcos,...
                                                     self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2),self.r_thetathetatheta_cyl(:,2),self.d_r_theta2_cyl_dRcos, self.d_r_thetatheta2_cyl_dRcos, self.d_r_thetathetatheta2_cyl_dRcos,...
                                                     self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3),self.r_thetathetatheta_cyl(:,3),self.d_r_theta3_cyl_dRcos, self.d_r_thetatheta3_cyl_dRcos, self.d_r_thetathetatheta3_cyl_dRcos,...
                                                     self.zeta, self.cR);

            self.d_torsion_dZsin = compute_torsion_dc(self.torsion_numerator,self.torsion_denominator, ...
                                                     self.r_theta_cyl(:,1), self.r_thetatheta_cyl(:,1),self.r_thetathetatheta_cyl(:,1),self.d_r_theta1_cyl_dZsin, self.d_r_thetatheta1_cyl_dZsin, self.d_r_thetathetatheta1_cyl_dZsin,...
                                                     self.r_theta_cyl(:,2), self.r_thetatheta_cyl(:,2),self.r_thetathetatheta_cyl(:,2),self.d_r_theta2_cyl_dZsin, self.d_r_thetatheta2_cyl_dZsin, self.d_r_thetathetatheta2_cyl_dZsin,...
                                                     self.r_theta_cyl(:,3), self.r_thetatheta_cyl(:,3),self.r_thetathetatheta_cyl(:,3),self.d_r_theta3_cyl_dZsin, self.d_r_thetatheta3_cyl_dZsin, self.d_r_thetathetatheta3_cyl_dZsin,...
                                                     self.zeta, self.sZ);
        
        

        
            self.tau_kappa2 = self.torsion./self.curvature.^2;
            self.d_tau_kappa2_dRcos = (self.curvature.^2 .* self.d_torsion_dRcos - 2 .* self.torsion .* self.curvature .* self.d_curvature_dRcos)./self.curvature.^4;
            self.d_tau_kappa2_dZsin = (self.curvature.^2 .* self.d_torsion_dZsin - 2 .* self.torsion .* self.curvature .* self.d_curvature_dZsin)./self.curvature.^4;

            self.d_abs_G0B0taukappa2_dRcos =  self.tau_kappa2.*self.d_abs_G_over_B0_dRcos+ self.abs_G_over_B0.*self.d_tau_kappa2_dRcos;
            self.d_abs_G0B0taukappa2_dZsin =  self.tau_kappa2.*self.d_abs_G_over_B0_dZsin+ self.abs_G_over_B0.*self.d_tau_kappa2_dZsin;
        
            % ------------------------------------------------------------------------------------
            % On axis field and gradient of field
            % ------------------------------------------------------------------------------------


            self.dBx_dRcos = self.B0 * self.dt1_dRcos;
            self.dBy_dRcos = self.B0 * self.dt2_dRcos;
            self.dBz_dRcos = self.B0 * self.dt3_dRcos;






            self.dBx_dZsin = self.B0 * self.dt1_dZsin;
            self.dBy_dZsin = self.B0 * self.dt2_dZsin;
            self.dBz_dZsin = self.B0 * self.dt3_dZsin;


            self.d_X1c_dRcos = -self.eta_bar *self.d_curvature_dRcos./ self.curvature.^2;
            self.d_X1c_dZsin = -self.eta_bar *self.d_curvature_dZsin./ self.curvature.^2;


            self.d_Y1s_dRcos =  self.d_curvature_dRcos./ self.eta_bar;
            self.d_Y1s_dZsin =  self.d_curvature_dZsin./ self.eta_bar;


            self.d_Y1c_dRcos =  self.d_curvature_dRcos.* self.sigma ./ self.eta_bar;
            self.d_Y1c_dZsin =  self.d_curvature_dZsin.* self.sigma ./ self.eta_bar;



            self.d2_X1c_dRcos = compute_d2_X1c_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dRcos,...
                                             self.s_theta, self.d_s_theta_dRcos,...
                                             self.curvature, self.d_curvature_dRcos, self.d_curvature_dtheta, self.d_curvature_dtheta_dRcos,...
                                             self.zeta, self.cR, self.eta_bar);
            self.d2_X1c_dZsin = compute_d2_X1c_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dZsin,...
                                             self.s_theta, self.d_s_theta_dZsin,...
                                             self.curvature, self.d_curvature_dZsin, self.d_curvature_dtheta, self.d_curvature_dtheta_dZsin,...
                                             self.zeta, self.sZ, self.eta_bar);


            self.d2_Y1s_dRcos = compute_d2_Y1s_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dRcos,...
                                             self.s_theta, self.d_s_theta_dRcos,...
                                             self.curvature,  self.d_curvature_dtheta, self.d_curvature_dtheta_dRcos,...
                                             self.zeta, self.cR, self.eta_bar);
            
            self.d2_Y1s_dZsin = compute_d2_Y1s_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dZsin,...
                                             self.s_theta, self.d_s_theta_dZsin,...
                                             self.curvature, self.d_curvature_dtheta, self.d_curvature_dtheta_dZsin,...
                                             self.zeta, self.sZ, self.eta_bar);
                             
            self.d2_Y1c_dRcos = compute_d2_Y1c_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dRcos,...
                                             self.s_theta, self.d_s_theta_dRcos,...
                                             self.curvature, self.d_curvature_dRcos, self.d_curvature_dtheta, self.d_curvature_dtheta_dRcos,...
                                             self.zeta, self.cR, self.eta_bar, self.D, self.sigma);
                                         

                                         
                                         
            self.d2_Y1c_dZsin = compute_d2_Y1c_dc(self.abs_G_over_B0, self.d_abs_G_over_B0_dZsin,...
                                             self.s_theta, self.d_s_theta_dZsin,...
                                             self.curvature, self.d_curvature_dZsin, self.d_curvature_dtheta, self.d_curvature_dtheta_dZsin,...
                                             self.zeta, self.sZ, self.eta_bar, self.D, self.sigma);
        
        
            self.d_X1c_detabar = 1./self.curvature;
            self.d_Y1s_detabar = self.svp * self.sg .* self.curvature .* (-1./self.eta_bar.^2);
            self.d_Y1c_detabar = self.svp * self.sg .* self.curvature .* (-1./self.eta_bar.^2) .* self.sigma;
        
            self.d2_X1c_detabar = self.d_X1c / self.eta_bar;
            self.d2_Y1s_detabar = self.eta_bar*self.d_Y1s * (-1./self.eta_bar.^2);
            self.d2_Y1c_detabar = self.eta_bar*self.d_Y1c * (-1./self.eta_bar.^2);
        
        
            [self.d_gradB1_detabar, self.d_gradB2_detabar, self.d_gradB3_detabar] = compute_gradB_detabar(self.B0, self.abs_G_over_B0,...
                                                                                                          self.iota, ...
                                                                                                          self.t_vec,self.n_vec,self.b_vec,...
                                                                                                          self.X1c, self.Y1s, self.Y1c, self.d_X1c, self.d_Y1s, self.d_Y1c,...
                                                                                                          self.d_X1c_detabar, self.d_Y1s_detabar, self.d_Y1c_detabar, ...
                                                                                                          self.d2_X1c_detabar, self.d2_Y1s_detabar, self.d2_Y1c_detabar);
        

                                     
        [self.d_B1x_dRcos, self.d_B1y_dRcos, self.d_B1z_dRcos,...
         self.d_B2x_dRcos, self.d_B2y_dRcos, self.d_B2z_dRcos,...
         self.d_B3x_dRcos, self.d_B3y_dRcos, self.d_B3z_dRcos] = compute_gradB_dc(self.abs_G_over_B0,self.d_abs_G_over_B0_dRcos,...
                                                                   self.torsion, self.curvature, ...
                                                                   self.d_torsion_dRcos, self.d_curvature_dRcos,...
                                                                   self.t_vec,self.n_vec,self.b_vec,...
                                                                   self.dt1_dRcos, self.dn1_dRcos,self.db1_dRcos,...
                                                                   self.dt2_dRcos, self.dn2_dRcos,self.db2_dRcos,...
                                                                   self.dt3_dRcos, self.dn3_dRcos,self.db3_dRcos,...
                                                                   self.X1c, self.Y1s, self.Y1c,self.d_X1c_dRcos, self.d_Y1s_dRcos, self.d_Y1c_dRcos,...
                                                                   self.d_X1c, self.d_Y1s, self.d_Y1c, self.d2_X1c_dRcos,self.d2_Y1s_dRcos,self.d2_Y1c_dRcos,...
                                                                   self.s_theta, self.d_s_theta_dRcos,...
                                                                   self.B0, self.iota, self.zeta, self.cR);
                                                                                                                     
                                                               

                                                               
                                                               
                                                               
        [self.d_B1x_dZsin, self.d_B1y_dZsin, self.d_B1z_dZsin,...
         self.d_B2x_dZsin, self.d_B2y_dZsin, self.d_B2z_dZsin,...
         self.d_B3x_dZsin, self.d_B3y_dZsin, self.d_B3z_dZsin] = compute_gradB_dc(self.abs_G_over_B0,self.d_abs_G_over_B0_dZsin,...
                                                                   self.torsion, self.curvature, ...
                                                                   self.d_torsion_dZsin, self.d_curvature_dZsin,...
                                                                   self.t_vec,self.n_vec,self.b_vec,...
                                                                   self.dt1_dZsin, self.dn1_dZsin,self.db1_dZsin,...
                                                                   self.dt2_dZsin, self.dn2_dZsin,self.db2_dZsin,...
                                                                   self.dt3_dZsin, self.dn3_dZsin,self.db3_dZsin,...
                                                                   self.X1c, self.Y1s, self.Y1c,self.d_X1c_dZsin, self.d_Y1s_dZsin, self.d_Y1c_dZsin,...
                                                                   self.d_X1c, self.d_Y1s, self.d_Y1c, self.d2_X1c_dZsin, self.d2_Y1s_dZsin, self.d2_Y1c_dZsin,...
                                                                   self.s_theta, self.d_s_theta_dZsin,...
                                                                   self.B0, self.iota, self.zeta, self.sZ);
        
        
        
        [self.d_B1x_ds, self.d_B1y_ds, self.d_B1z_ds,...
         self.d_B2x_ds, self.d_B2y_ds, self.d_B2z_ds,...
         self.d_B3x_ds, self.d_B3y_ds, self.d_B3z_ds] = compute_gradB_ds(self.B0, self.abs_G_over_B0,...
                                                                         self.zeta, self.sigma, self.iota,...
                                                                         self.n_vec,self.b_vec,...
                                                                         self.X1c, self.Y1s, self.Y1c,...
                                                                         self.d_Y1c_ds,self.d_Y1s, self.d2_Y1c_ds);
                                                          
        [self.d_gradB1_diota, self.d_gradB2_diota, self.d_gradB3_diota] = compute_gradB_diota(self.B0, self.abs_G_over_B0,...
                                                                                  self.n_vec,self.b_vec, ...
                                                                                  self.X1c, self.Y1s, self.Y1c);
        
        

                                                      
                                                                              
                                                                              
                                                      
                                                                              
        end
        


        
        function [] = compute_QS_axis_data(self)
            %       (t,n,b) and its derivatives in Cartesian coordinates
            t_cyl = self.r_theta_cyl./self.s_theta;
            n_cyl = self.kappa./self.curvature;
            self.t_vec = [cos(self.zeta).*t_cyl(:,1)-sin(self.zeta).*t_cyl(:,2), sin(self.zeta).*t_cyl(:,1)+cos(self.zeta).*t_cyl(:,2), t_cyl(:,3)];
            self.n_vec = [cos(self.zeta).*n_cyl(:,1)-sin(self.zeta).*n_cyl(:,2), sin(self.zeta).*n_cyl(:,1)+cos(self.zeta).*n_cyl(:,2), n_cyl(:,3)];
            self.b_vec = cross(self.t_vec,self.n_vec);
        
            
            
            % ------------------------------------------------------------------------------------
            % On axis field and gradient of field
            % ------------------------------------------------------------------------------------

            self.X1c = self.eta_bar ./ self.curvature;
            self.Y1s = self.curvature ./ self.eta_bar;
            self.Y1c = self.curvature .* self.sigma / self.eta_bar;

            self.d_X1c = (self.abs_G_over_B0./self.s_theta).*self.eta_bar .* (-self.d_curvature_dtheta./self.curvature.^2);
            self.d_Y1s = (self.abs_G_over_B0./self.s_theta)./self.eta_bar .* self.d_curvature_dtheta;
            self.d_Y1c = (self.abs_G_over_B0./self.s_theta)./self.eta_bar .* (self.d_curvature_dtheta.*self.sigma + (self.D*self.sigma).*self.curvature);

            self.d_Y1c_ds = diag(self.curvature ./ self.eta_bar);
            self.d2_Y1c_ds = (self.abs_G_over_B0./self.s_theta)./self.eta_bar .* (diag(self.d_curvature_dtheta) + self.D .* self.curvature);

            
            self.B_QS = self.B0*self.t_vec;
            
            self.gradBx_QS = self.B0./self.abs_G_over_B0.*((self.abs_G_over_B0 .* self.curvature .* self.t_vec(:,1) + (self.d_X1c.*self.Y1s+self.iota*self.X1c.*self.Y1c).*self.n_vec(:,1) + ...
                                         (self.d_Y1c.*self.Y1s-self.d_Y1s.*self.Y1c+self.abs_G_over_B0.*self.torsion+self.iota*(self.Y1s.^2 + self.Y1c.^2) ).*self.b_vec(:,1) ) .* self.n_vec + ...
                                         ( (-self.abs_G_over_B0.*self.torsion-self.iota*self.X1c.^2).*self.n_vec(:,1)+(self.X1c.*self.d_Y1s-self.iota*self.X1c.*self.Y1c).*self.b_vec(:,1) ) .* self.b_vec ) +...
                                         self.curvature.*self.B0.*self.n_vec(:,1).*self.t_vec;
            self.gradBy_QS = self.B0./self.abs_G_over_B0.*((self.abs_G_over_B0 .* self.curvature .* self.t_vec(:,2) + (self.d_X1c.*self.Y1s+self.iota*self.X1c.*self.Y1c).*self.n_vec(:,2) + ...
                                         (self.d_Y1c.*self.Y1s-self.d_Y1s.*self.Y1c+self.abs_G_over_B0.*self.torsion+self.iota*(self.Y1s.^2 + self.Y1c.^2) ).*self.b_vec(:,2) ) .* self.n_vec + ...
                                         ( (-self.abs_G_over_B0.*self.torsion-self.iota*self.X1c.^2).*self.n_vec(:,2)+(self.X1c.*self.d_Y1s-self.iota*self.X1c.*self.Y1c).*self.b_vec(:,2) ) .* self.b_vec ) +...
                                         self.curvature.*self.B0.*self.n_vec(:,2).*self.t_vec;
            self.gradBz_QS = self.B0./self.abs_G_over_B0.*((self.abs_G_over_B0 .* self.curvature .* self.t_vec(:,3) + (self.d_X1c.*self.Y1s+self.iota*self.X1c.*self.Y1c).*self.n_vec(:,3) + ...
                                         (self.d_Y1c.*self.Y1s-self.d_Y1s.*self.Y1c+self.abs_G_over_B0.*self.torsion+self.iota*(self.Y1s.^2 + self.Y1c.^2) ).*self.b_vec(:,3) ) .* self.n_vec + ...
                                         ( (-self.abs_G_over_B0.*self.torsion-self.iota*self.X1c.^2).*self.n_vec(:,3)+(self.X1c.*self.d_Y1s-self.iota*self.X1c.*self.Y1c).*self.b_vec(:,3) ) .* self.b_vec ) +...
                                         self.curvature.*self.B0.*self.n_vec(:,3).*self.t_vec; 


        end
        
        
        function [] = updateMagneticAxisFrameData(self)
            self.RA = zeros(self.nzeta,1);
            self.ZA = zeros(self.nzeta,1);
            self.RAp = zeros(self.nzeta,1);
            self.ZAp = zeros(self.nzeta,1);
            self.RApp = zeros(self.nzeta,1);
            self.ZApp = zeros(self.nzeta,1);
            self.RAppp = zeros(self.nzeta,1);
            self.ZAppp = zeros(self.nzeta,1);
            
            % compute frame data for magnetic axis
            for imn = 1:self.top
                n = self.nfp*(imn-1);
                angle = n * self.zeta;
                sinangle = sin(angle);
                cosangle = cos(angle);
                self.RA = self.RA + self.cR(imn)*cosangle ;
                self.ZA = self.ZA + self.sZ(imn)*sinangle;
                self.RAp = self.RAp - n*self.cR(imn)*sinangle ;
                self.ZAp = self.ZAp + n*self.sZ(imn)*cosangle;
                self.RApp = self.RApp - n*n*self.cR(imn)*cosangle ;
                self.ZApp = self.ZApp - n*n*self.sZ(imn)*sinangle;
                self.RAppp = self.RAppp + n*n*n*self.cR(imn)*sinangle ;
                self.ZAppp = self.ZAppp - n*n*n*self.sZ(imn)*cosangle;
            end
            
            
            
            dR_dPhi =  self.dCOS*self.cR;

            self.dX_dPhi = dR_dPhi .* cos(self.zeta) -  self.RA .* sin(self.zeta);
            self.dY_dPhi = dR_dPhi .* sin(self.zeta) +  self.RA .* cos(self.zeta);
            self.dZ_dPhi = self.dSIN*self.sZ ;
            
            self.Ta_mag = sqrt(self.dX_dPhi.^2 + self.dY_dPhi.^2 + self.dZ_dPhi.^2);
            
            
            
            self.s_theta = sqrt(self.RA .* self.RA + self.RAp .* self.RAp + self.ZAp .* self.ZAp);
            self.s_thetatheta = (self.RA .* self.RAp + self.RAp .* self.RApp + self.ZAp .* self.ZApp) ./ self.s_theta;
            self.s_thetathetatheta = (self.RA .* self.RApp + self.RAp .^ 2 + self.RAp .* self.RAppp + self.RApp .^ 2 + self.ZAp .* self.ZAppp + self.ZApp .^ 2) ./ self.s_theta ...
                              - (self.RA .* self.RAp + self.RAp .* self.RApp + self.ZAp .* self.ZApp) ./ self.s_theta .^ 2 .* self.s_thetatheta;
            
             
            self.r_theta_cyl = [self.RAp, self.RA, self.ZAp]; 
            self.r_thetatheta_cyl = [self.RApp-self.RA, 2*self.RAp, self.ZApp];  
            self.r_thetathetatheta_cyl = [self.RAppp-3*self.RAp, 3*self.RApp-self.RA, self.ZAppp];
        
            
            self.r_theta_cyl_tc = [self.RApp, self.RAp, self.ZApp];
            self.r_thetatheta_cyl_tc = [self.RAppp-self.RAp, 2*self.RApp, self.ZAppp]; 

            self.kappa = [...
                (-self.r_theta_cyl(:,1) .* self.s_thetatheta ./ self.s_theta + self.r_thetatheta_cyl(:,1)) ./ (self.s_theta .* self.s_theta), ...
                (-self.r_theta_cyl(:,2) .* self.s_thetatheta ./ self.s_theta + self.r_thetatheta_cyl(:,2)) ./ (self.s_theta .* self.s_theta), ...
                (-self.r_theta_cyl(:,3) .* self.s_thetatheta ./ self.s_theta + self.r_thetatheta_cyl(:,3)) ./ (self.s_theta .* self.s_theta)];

            self.curvature = sqrt(self.kappa(:,1) .* self.kappa(:,1) + self.kappa(:,2) .* self.kappa(:,2) + self.kappa(:,3) .* self.kappa(:,3));

            self.d_kappa_dtheta1 = (self.r_thetatheta_cyl_tc(:,1) .* self.s_theta .^ 2 + ((-self.r_theta_cyl_tc(:,1) - 2 * self.r_thetatheta_cyl(:,1)) .* self.s_thetatheta - self.r_theta_cyl(:,1) .* self.s_thetathetatheta) .* self.s_theta + 3 * self.r_theta_cyl(:,1) .* self.s_thetatheta .^ 2) ./ self.s_theta .^ 4;
            self.d_kappa_dtheta2 = (self.r_thetatheta_cyl_tc(:,2) .* self.s_theta .^ 2 + ((-self.r_theta_cyl_tc(:,2) - 2 * self.r_thetatheta_cyl(:,2)) .* self.s_thetatheta - self.r_theta_cyl(:,2) .* self.s_thetathetatheta) .* self.s_theta + 3 * self.r_theta_cyl(:,2) .* self.s_thetatheta .^ 2) ./ self.s_theta .^ 4;
            self.d_kappa_dtheta3 = (self.r_thetatheta_cyl_tc(:,3) .* self.s_theta .^ 2 + ((-self.r_theta_cyl_tc(:,3) - 2 * self.r_thetatheta_cyl(:,3)) .* self.s_thetatheta - self.r_theta_cyl(:,3) .* self.s_thetathetatheta) .* self.s_theta + 3 * self.r_theta_cyl(:,3) .* self.s_thetatheta .^ 2) ./ self.s_theta .^ 4;
        
            self.d_kappa_dtheta = [self.d_kappa_dtheta1, self.d_kappa_dtheta2, self.d_kappa_dtheta3];
            self.d_curvature_dtheta = (self.kappa(:,1) .* self.d_kappa_dtheta(:,1) + self.kappa(:,2) .* self.d_kappa_dtheta(:,2) + self.kappa(:,3) .* self.d_kappa_dtheta(:,3))./self.curvature;
        
            

            
            % This script uses the same sign convention for torsion as the J Plasma Phys papers, wikipedia
            % and mathworld.wolfram.com/Torsion.html, opposite to Garren & Boozer's sign convention.
            self.torsion_numerator = (0 ...
                + self.r_theta_cyl(:,1) .* (self.r_thetatheta_cyl(:,2) .* self.r_thetathetatheta_cyl(:,3) - self.r_thetatheta_cyl(:,3) .* self.r_thetathetatheta_cyl(:,2)) ...
                + self.r_theta_cyl(:,2) .* (self.r_thetatheta_cyl(:,3) .* self.r_thetathetatheta_cyl(:,1) - self.r_thetatheta_cyl(:,1) .* self.r_thetathetatheta_cyl(:,3)) ...
                + self.r_theta_cyl(:,3) .* (self.r_thetatheta_cyl(:,1) .* self.r_thetathetatheta_cyl(:,2) - self.r_thetatheta_cyl(:,2) .* self.r_thetathetatheta_cyl(:,1)));

            self.torsion_denominator = 0 ...
                + (self.r_theta_cyl(:,2) .* self.r_thetatheta_cyl(:,3) - self.r_theta_cyl(:,3) .* self.r_thetatheta_cyl(:,2)) .^ 2 ...
                + (self.r_theta_cyl(:,3) .* self.r_thetatheta_cyl(:,1) - self.r_theta_cyl(:,1) .* self.r_thetatheta_cyl(:,3)) .^ 2 ...
                + (self.r_theta_cyl(:,1) .* self.r_thetatheta_cyl(:,2) - self.r_theta_cyl(:,2) .* self.r_thetatheta_cyl(:,1)) .^ 2;

            self.torsion = self.torsion_numerator ./ self.torsion_denominator;
            self.Dtilde = diag(mean(self.s_theta) ./ (self.B0* self.s_theta)) * self.D;
            
            self.I2_over_B0          = 0.;
            self.abs_G_over_B0       = mean(self.s_theta)/self.B0;  
        end
        
        function min_residual = determineSigmaAndIota(self)
            min_residual = 0;
            
            self.updateMagneticAxisFrameData();
            state = [self.sigma; self.iota];
            options=optimoptions('fsolve','Display','off', 'SpecifyObjectiveGradient',true, 'OptimalityTolerance', 1e-12);
%             options=optimoptions('fsolve','Display','iter-detailed', 'SpecifyObjectiveGradient',true, 'OptimalityTolerance', 1e-12);  
            

            [state, min_residual] = fsolve(@(x)compute_FJ(x, self.eta_bar, self.curvature, self.torsion, self.Dtilde, self.I2_over_B0, self.abs_G_over_B0), state, options);
            self.sigma = state(1:end-1);
            self.iota = state(end);
            
            
        end
        
        
        % we need B and gradB on axis generated by the coils.
        function compute_coil_axis_data(self, coilData, adjoint)
            
            self.B_coils = zeros(self.nzeta,3);
            self.gradBx_coils = zeros(self.nzeta,3);
            self.gradBy_coils = zeros(self.nzeta,3);
            self.gradBz_coils = zeros(self.nzeta,3);

            if adjoint
                fun = @(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9)biotsavartwithpartials(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);

                for k = 1: self.nzeta
                    [self.B_coils(k,1), self.B_coils(k,2), self.B_coils(k,3),...
                        self.gradBx_coils(k,1), self.gradBx_coils(k,2), self.gradBx_coils(k,3),...
                        self.gradBy_coils(k,1), self.gradBy_coils(k,2), self.gradBy_coils(k,3),...
                        self.gradBz_coils(k,1), self.gradBz_coils(k,2), self.gradBz_coils(k,3),...
                        self.dBx_coils_dc(k,:), self.dBy_coils_dc(k,:), self.dBz_coils_dc(k,:),...
                        self.dBx_coils_dI(k,:), self.dBy_coils_dI(k,:), self.dBz_coils_dI(k,:),...
                        self.dB1x_dI(k,:), self.dB1y_dI(k,:), self.dB1z_dI(k,:),...
                        self.dB2x_dI(k,:), self.dB2y_dI(k,:), self.dB2z_dI(k,:),...
                        self.dB3x_dI(k,:), self.dB3y_dI(k,:), self.dB3z_dI(k,:),...
                        self.dBx_coils_dRcos(k,:), self.dBy_coils_dRcos(k,:), self.dBz_coils_dRcos(k,:),...
                        self.dBx_coils_dZsin(k,:), self.dBy_coils_dZsin(k,:), self.dBz_coils_dZsin(k,:),...
                        self.dB1x_dRcos(k,:),self.dB1y_dRcos(k,:),self.dB1z_dRcos(k,:),...
                        self.dB2x_dRcos(k,:),self.dB2y_dRcos(k,:),self.dB2z_dRcos(k,:),...
                        self.dB3x_dRcos(k,:),self.dB3y_dRcos(k,:),self.dB3z_dRcos(k,:),...
                        self.dB1x_dZsin(k,:),self.dB1y_dZsin(k,:),self.dB1z_dZsin(k,:),...
                        self.dB2x_dZsin(k,:),self.dB2y_dZsin(k,:),self.dB2z_dZsin(k,:),...
                        self.dB3x_dZsin(k,:),self.dB3y_dZsin(k,:),self.dB3z_dZsin(k,:),...
                        self.B1x_coils_dc(k,:), self.B1y_coils_dc(k,:), self.B1z_coils_dc(k,:), ...
                        self.B2x_coils_dc(k,:), self.B2y_coils_dc(k,:), self.B2z_coils_dc(k,:), ...
                        self.B3x_coils_dc(k,:), self.B3y_coils_dc(k,:), self.B3z_coils_dc(k,:)] = fun(self.RA(k)*cos(self.zeta(k)), self.RA(k)*sin(self.zeta(k)),self.ZA(k),coilData,...
                            self.dRA_dcos(k,:)*cos(self.zeta(k)), self.dRA_dcos(k,:)*sin(self.zeta(k)), zeros(1,self.top),...
                            zeros(1,self.top), zeros(1,self.top), self.dZA_dsin(k,:));
                end
                
            else
                for k = 1: self.nzeta
                    [self.B_coils(k,1), self.B_coils(k,2), self.B_coils(k,3),...
                     self.gradBx_coils(k,1), self.gradBx_coils(k,2), self.gradBx_coils(k,3),...
                     self.gradBy_coils(k,1), self.gradBy_coils(k,2), self.gradBy_coils(k,3),...
                     self.gradBz_coils(k,1), self.gradBz_coils(k,2), self.gradBz_coils(k,3)] = biotsavartAndGrad(self.RA(k)*cos(self.zeta(k)), self.RA(k)*sin(self.zeta(k)),self.ZA(k),coilData);
                end
            end
        end
        
        function B_diff = computeBdifference(self)
            dB_diff  = self.B_QS-self.B_coils;
            dB_diff2 = dB_diff(:,1).^2+dB_diff(:,2).^2+dB_diff(:,3).^2;
            B_diff   = 0.5 * (2*pi/length(self.zeta))*sum( self.Ta_mag.*dB_diff2 );
        end
        
        function [dBdiff_dcc] = compute_dBdiff_dcc(self)
            dB_diff = self.B_QS-self.B_coils;
            temp = -(dB_diff(:,1) .* self.dBx_coils_dc + dB_diff(:,2) .* self.dBy_coils_dc + dB_diff(:,3) .* self.dBz_coils_dc).*self.Ta_mag;
            dBdiff_dcc = (2*pi/length(self.zeta))*sum(temp,1);
        end
        
        function [dBdiff_dca] = compute_dBdiff_dca(self)

            dB_diff = self.B_QS-self.B_coils;
            
            temp = (dB_diff(:,1) .* (self.dBx_dRcos-self.dBx_coils_dRcos) + dB_diff(:,2) .* (self.dBy_dRcos-self.dBy_coils_dRcos) + dB_diff(:,3) .* (self.dBz_dRcos-self.dBz_coils_dRcos) ).*self.Ta_mag...
                  +(dB_diff(:,1).^2+dB_diff(:,2).^2+dB_diff(:,3).^2).*self.d_Ta_mag_dRcos*0.5;
            dBdiff_dRcos = (2*pi/length(self.zeta))*sum(temp,1);
            

            
            
            temp = (dB_diff(:,1) .* (self.dBx_dZsin-self.dBx_coils_dZsin) + dB_diff(:,2) .* (self.dBy_dZsin-self.dBy_coils_dZsin) + dB_diff(:,3) .* (self.dBz_dZsin-self.dBz_coils_dZsin) ).*self.Ta_mag...
                  +(dB_diff(:,1).^2+dB_diff(:,2).^2+dB_diff(:,3).^2).*self.d_Ta_mag_dZsin*0.5;
            dBdiff_dZsin = (2*pi/length(self.zeta))*sum(temp,1);
            
            dBdiff_dca = [dBdiff_dRcos dBdiff_dZsin];
        end
        
        
        function gradB_diff = computeGradBdifference(self)
            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
            
            gradB_diff2 =  gradBx_diff(:,1).^2+gradBx_diff(:,2).^2+gradBx_diff(:,3).^2 ...
                         + gradBy_diff(:,1).^2+gradBy_diff(:,2).^2+gradBy_diff(:,3).^2 ...
                         + gradBz_diff(:,1).^2+gradBz_diff(:,2).^2+gradBz_diff(:,3).^2;
            gradB_diff = 0.5 * (2*pi/length(self.zeta))*sum( self.Ta_mag.*gradB_diff2 );
        end
        
        function [dBdiff_dcc] = compute_dGradBdiff_dcc(self)
            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
            
            temp = -(gradBx_diff(:,1) .* self.B1x_coils_dc + gradBx_diff(:,2) .* self.B1y_coils_dc + gradBx_diff(:,3) .* self.B1z_coils_dc + ...
                     gradBy_diff(:,1) .* self.B2x_coils_dc + gradBy_diff(:,2) .* self.B2y_coils_dc + gradBy_diff(:,3) .* self.B2z_coils_dc + ...
                     gradBz_diff(:,1) .* self.B3x_coils_dc + gradBz_diff(:,2) .* self.B3y_coils_dc + gradBz_diff(:,3) .* self.B3z_coils_dc).*self.Ta_mag;
            dBdiff_dcc = (2*pi/self.nzeta)*sum(temp,1);
        end
        
        function [dBdiff_dca] = compute_dGradBdiff_dca(self)

            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
 
            temp =  (gradBx_diff(:,1) .* (self.d_B1x_dRcos-self.dB1x_dRcos) + gradBx_diff(:,2) .* (self.d_B1y_dRcos-self.dB1y_dRcos) + gradBx_diff(:,3) .* (self.d_B1z_dRcos-self.dB1z_dRcos) ...
                    +gradBy_diff(:,1) .* (self.d_B2x_dRcos-self.dB2x_dRcos) + gradBy_diff(:,2) .* (self.d_B2y_dRcos-self.dB2y_dRcos) + gradBy_diff(:,3) .* (self.d_B2z_dRcos-self.dB2z_dRcos) ...
                    +gradBz_diff(:,1) .* (self.d_B3x_dRcos-self.dB3x_dRcos) + gradBz_diff(:,2) .* (self.d_B3y_dRcos-self.dB3y_dRcos) + gradBz_diff(:,3) .* (self.d_B3z_dRcos-self.dB3z_dRcos) ) .* self.Ta_mag...
                    +(gradBx_diff(:,1).^2+gradBx_diff(:,2).^2+gradBx_diff(:,3).^2 ...
                     +gradBy_diff(:,1).^2+gradBy_diff(:,2).^2+gradBy_diff(:,3).^2 ...
                     +gradBz_diff(:,1).^2+gradBz_diff(:,2).^2+gradBz_diff(:,3).^2).*self.d_Ta_mag_dRcos*0.5;

            dgradBdiff_dRcos = (2*pi/length(self.zeta))*sum(temp,1);
            

            temp =  (gradBx_diff(:,1) .* (self.d_B1x_dZsin-self.dB1x_dZsin) + gradBx_diff(:,2) .* (self.d_B1y_dZsin-self.dB1y_dZsin) + gradBx_diff(:,3) .* (self.d_B1z_dZsin-self.dB1z_dZsin) ...
                    +gradBy_diff(:,1) .* (self.d_B2x_dZsin-self.dB2x_dZsin) + gradBy_diff(:,2) .* (self.d_B2y_dZsin-self.dB2y_dZsin) + gradBy_diff(:,3) .* (self.d_B2z_dZsin-self.dB2z_dZsin) ...
                    +gradBz_diff(:,1) .* (self.d_B3x_dZsin-self.dB3x_dZsin) + gradBz_diff(:,2) .* (self.d_B3y_dZsin-self.dB3y_dZsin) + gradBz_diff(:,3) .* (self.d_B3z_dZsin-self.dB3z_dZsin) ) .* self.Ta_mag...
                    +(gradBx_diff(:,1).^2+gradBx_diff(:,2).^2+gradBx_diff(:,3).^2 ...
                     +gradBy_diff(:,1).^2+gradBy_diff(:,2).^2+gradBy_diff(:,3).^2 ...
                     +gradBz_diff(:,1).^2+gradBz_diff(:,2).^2+gradBz_diff(:,3).^2).*self.d_Ta_mag_dZsin*0.5;

            dgradBdiff_dZsin = (2*pi/length(self.zeta))*sum(temp,1);
            dBdiff_dca = [dgradBdiff_dRcos dgradBdiff_dZsin];
        end
        
        function dgradBdiff_ds = compute_gradBdiff_ds(self)
%             dgradBdiff_ds = zeros(1,self.nzeta + 1);


            temp_sigma = zeros(1,self.nzeta);
            temp_iota = 0.;
            
            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
            
            for np = 1:self.nzeta
                temp_sigma =  temp_sigma + ...
                        (gradBx_diff(np,1) * self.d_B1x_ds(np,:) + gradBx_diff(np,2) * self.d_B1y_ds(np,:) + gradBx_diff(np,3) * self.d_B1z_ds(np,:) ...
                        +gradBy_diff(np,1) * self.d_B2x_ds(np,:) + gradBy_diff(np,2) * self.d_B2y_ds(np,:) + gradBy_diff(np,3) * self.d_B2z_ds(np,:) ...
                        +gradBz_diff(np,1) * self.d_B3x_ds(np,:) + gradBz_diff(np,2) * self.d_B3y_ds(np,:) + gradBz_diff(np,3) * self.d_B3z_ds(np,:) ) * self.Ta_mag(np);
            end
            temp_sigma = (2*pi/self.nzeta)*temp_sigma;  
            
            for np = 1:self.nzeta
                temp_iota =  temp_iota + ...
	                        (gradBx_diff(np,1) * self.d_gradB1_diota(np,1) + gradBx_diff(np,2) * self.d_gradB1_diota(np,2) + gradBx_diff(np,3) * self.d_gradB1_diota(np,3) ...
	                        +gradBy_diff(np,1) * self.d_gradB2_diota(np,1) + gradBy_diff(np,2) * self.d_gradB2_diota(np,2) + gradBy_diff(np,3) * self.d_gradB2_diota(np,3) ...
	                        +gradBz_diff(np,1) * self.d_gradB3_diota(np,1) + gradBz_diff(np,2) * self.d_gradB3_diota(np,2) + gradBz_diff(np,3) * self.d_gradB3_diota(np,3) ) * self.Ta_mag(np);
            end
            temp_iota = (2*pi/self.nzeta)*temp_iota;  
            
            
            dgradBdiff_ds = [temp_sigma temp_iota];
                     
                     
%             dgradBdiff_ds = 0.5 * (2*pi/length(self.zeta))*sum( self.Ta_mag.*gradB_diff2 );
        end
        
       function [len] = get_alength(self)
            len = (2*pi/length(self.zeta))*sum( self.Ta_mag );
        end 
        
        function [dalength_dca] = compute_dalength_dca(self)
            dalength_dRcos = (2*pi/length(self.zeta))*sum( self.d_Ta_mag_dRcos ,1);
            dalength_dZsin = (2*pi/length(self.zeta))*sum( self.d_Ta_mag_dZsin ,1);
            
            dalength_dca = [dalength_dRcos dalength_dZsin];
        end
        

         function [dBdiff_dI] = compute_dBdiff_dI(self)
            dB_diff = self.B_QS-self.B_coils;
            temp = -(dB_diff(:,1) .* self.dBx_coils_dI + dB_diff(:,2) .* self.dBy_coils_dI + dB_diff(:,3) .* self.dBz_coils_dI).*self.Ta_mag;
            dBdiff_dI = (2*pi/length(self.zeta))*sum(temp,1);
         end
        
         function [dGradBdiff_dI] = compute_dGradBdiff_dI(self)

            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
            
            temp = -(gradBx_diff(:,1) .* self.dB1x_dI + gradBx_diff(:,2) .* self.dB1y_dI + gradBx_diff(:,3) .* self.dB1z_dI + ...
                     gradBy_diff(:,1) .* self.dB2x_dI + gradBy_diff(:,2) .* self.dB2y_dI + gradBy_diff(:,3) .* self.dB2z_dI + ...
                     gradBz_diff(:,1) .* self.dB3x_dI + gradBz_diff(:,2) .* self.dB3y_dI + gradBz_diff(:,3) .* self.dB3z_dI).*self.Ta_mag;
            dGradBdiff_dI = (2*pi/self.nzeta)*sum(temp,1);

         end



         function [dGradBdiff_detabar] = compute_dGradBdiff_detabar(self)
             
            gradBx_diff = self.gradBx_QS-self.gradBx_coils;
            gradBy_diff = self.gradBy_QS-self.gradBy_coils;
            gradBz_diff = self.gradBz_QS-self.gradBz_coils;
            
            temp =  (gradBx_diff(:,1) .* self.d_gradB1_detabar(:,1) + gradBx_diff(:,2) .* self.d_gradB1_detabar(:,2) + gradBx_diff(:,3) .* self.d_gradB1_detabar(:,3) + ...
                     gradBy_diff(:,1) .* self.d_gradB2_detabar(:,1) + gradBy_diff(:,2) .* self.d_gradB2_detabar(:,2) + gradBy_diff(:,3) .* self.d_gradB2_detabar(:,3) + ...
                     gradBz_diff(:,1) .* self.d_gradB3_detabar(:,1) + gradBz_diff(:,2) .* self.d_gradB3_detabar(:,2) + gradBz_diff(:,3) .* self.d_gradB3_detabar(:,3)).*self.Ta_mag;
            dGradBdiff_detabar = (2*pi/self.nzeta)*sum(temp,1);
             
         end




        function matrix = compute_fca(self)  
            matrix = zeros(size(self.sigma,1)+1, numel(self.cR)+numel(self.sZ));
            offset = 0;
            for mode = 1:numel(self.cR)
                temp = (self.d_abs_G_over_B0_dRcos(:,mode).*self.s_theta - self.abs_G_over_B0.*self.d_s_theta_dRcos(:,mode))./self.s_theta.^2;
                matrix(:,offset + mode) = [diag(temp)*self.D*self.sigma - 4 * self.eta_bar.^4.*self.iota.*self.d_curvature_dRcos(:,mode)./self.curvature.^5  + 2 * self.eta_bar.^2 .* self.d_abs_G0B0taukappa2_dRcos(:,mode);0];
            end
            offset = offset + numel(self.cR);
%             for mode = 1:numel(self.sR)
%                 temp = (self.d_abs_G_over_B0_dRsin(:,mode).*self.s_theta - self.abs_G_over_B0.*self.d_s_theta_dRsin(:,mode))./self.s_theta.^2;
%                 matrix(:,offset + mode) = [diag(temp)*self.D*self.sigma - 4 * self.eta_bar.^4.*self.iota.*self.d_curvature_dRsin(:,mode)./self.curvature.^5 + 2 * self.eta_bar.^2 .* self.d_abs_G0B0taukappa2_dRsin(:,mode);0];
%             end
%             offset = offset + numel(self.sR);
%             for mode = 1:numel(self.cZ)
%                 temp = (self.d_abs_G_over_B0_dZcos(:,mode).*self.s_theta - self.abs_G_over_B0.*self.d_s_theta_dZcos(:,mode))./self.s_theta.^2;
%                 matrix(:,offset + mode) = [diag(temp)*self.D*self.sigma - 4 * self.eta_bar.^4.*self.iota.*self.d_curvature_dZcos(:,mode)./self.curvature.^5 + 2 * self.eta_bar.^2 .* self.d_abs_G0B0taukappa2_dZcos(:,mode);0];
%             end
%             offset = offset + numel(self.cZ);
            for mode = 1:numel(self.sZ)
                temp = (self.d_abs_G_over_B0_dZsin(:,mode).*self.s_theta - self.abs_G_over_B0.*self.d_s_theta_dZsin(:,mode))./self.s_theta.^2;
                matrix(:,offset + mode) = [diag(temp)*self.D*self.sigma  - 4 * self.eta_bar.^4.*self.iota.*self.d_curvature_dZsin(:,mode)./self.curvature.^5 + 2 * self.eta_bar.^2 .* self.d_abs_G0B0taukappa2_dZsin(:,mode);0];
            end
        end
        
        
        function [J] = compute_fs(self)
            J = zeros( size(self.sigma,1)+1 );
            J(1:end-1,end) =  (self.eta_bar*self.eta_bar*self.eta_bar*self.eta_bar) ./ (self.curvature .* self.curvature .* self.curvature .* self.curvature) + 1 + self.sigma .* self.sigma;
            J(1:end-1,1:end-1) = self.Dtilde + diag(self.iota * 2 * self.sigma);
            J(end, 1) = 1;
        end
        
        function [fetabar] = compute_fetabar(self)
            temp = 4 * self.iota * self.eta_bar.^3 ./ self.curvature.^4 ...
                 + 4 * self.abs_G_over_B0 .* self.eta_bar .* self.tau_kappa2;
            fetabar = [temp ; 0];
        end
        
        
     


        function [iota_ds] = compute_iota_ds(self)
            
            iota_ds = zeros(1,size(self.sigma,1)+1);
            iota_ds(end) = 1;
            
        end
        
        function [] = plot_magnetic_axis(self)
            self.updateMagneticAxisFrameData();
            
            x = self.RA.*cos(self.zeta);
            y = self.RA.*sin(self.zeta);
            z = self.ZA;
            for t = 1:self.nfp-1
                R = [cos(t * 2 * pi / self.nfp), - sin(t * 2 * pi / self.nfp) 0.; sin(t * 2 * pi / self.nfp), cos(t * 2 * pi / self.nfp), 0; 0 0 1];
                p = R * [x(1:self.nzeta) y(1:self.nzeta) z(1:self.nzeta)]';
                
                
                x = [x ; p(1,:)'] ;
                y = [y ; p(2,:)'] ;
                z = [z ; p(3,:)'] ;
            end
            x = [x ; x(1)] ;
            y = [y ; y(1)] ;
            z = [z ; z(1)];
            plot3(x, y, z, '-k'); hold on;
            
        end
        
        
        
    end


end

function [d_stt_dc] = compute_d_stt_dRc(R0, R0p, R0pp, R0ppp, Z0, Z0p, Z0pp, Z0ppp,...
                                        dR0_dc,dR0p_dc,dR0pp_dc,dR0ppp_dc,...
                                        s_theta, s_thetatheta, d_s_theta_dc, d_s_thetatheta_dc,...
                                        phi, R0c)
    d_stt_dc = zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d_stt_dc(:,mode)=(R0 .* dR0pp_dc(:,mode) + 2 .* R0p .* dR0p_dc(:,mode) + R0p .* dR0ppp_dc(:,mode) + dR0_dc(:,mode) .* R0pp + 2 .* R0pp .* dR0pp_dc(:,mode) + dR0p_dc(:,mode) .* R0ppp) ./ s_theta - (R0 .* R0pp + R0p .^ 2 + R0p .* R0ppp + R0pp .^ 2 + Z0p .* Z0ppp + Z0pp .^ 2) ./ s_theta .^ 2 .* d_s_theta_dc(:,mode) - (R0 .* dR0p_dc(:,mode) + dR0_dc(:,mode) .* R0p + R0p .* dR0pp_dc(:,mode) + dR0p_dc(:,mode) .* R0pp) ./ s_theta .^ 2 .* s_thetatheta + 2 .* (R0 .* R0p + R0p .* R0pp + Z0p .* Z0pp) ./ s_theta .^ 3 .* s_thetatheta .* d_s_theta_dc(:,mode) - (R0 .* R0p + R0p .* R0pp + Z0p .* Z0pp) ./ s_theta .^ 2 .* d_s_thetatheta_dc(:,mode);
    end
end
function [d_stt_dc] = compute_d_stt_dZc(R0, R0p, R0pp, R0ppp, Z0, Z0p, Z0pp, Z0ppp,...
                                        dZ0_dc, dZ0p_dc, dZ0pp_dc, dZ0ppp_dc,...
                                        s_theta, s_thetatheta, d_s_theta_dc, d_s_thetatheta_dc,...
                                        phi, R0c)
    d_stt_dc = zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d_stt_dc(:,mode)=(Z0p .* dZ0ppp_dc(:,mode) + 2 .* Z0pp .* dZ0pp_dc(:,mode) + Z0ppp .* dZ0p_dc(:,mode)) ./ s_theta - (R0 .* R0pp + R0p .^ 2 + R0p .* R0ppp + R0pp .^ 2 + Z0p .* Z0ppp + Z0pp .^ 2) ./ s_theta .^ 2 .* d_s_theta_dc(:,mode) - (Z0p .* dZ0pp_dc(:,mode) + Z0pp .* dZ0p_dc(:,mode)) ./ s_theta .^ 2 .* s_thetatheta + 2 .* (R0 .* R0p + R0p .* R0pp + Z0p .* Z0pp) ./ s_theta .^ 3 .* s_thetatheta .* d_s_theta_dc(:,mode) - (R0 .* R0p + R0p .* R0pp + Z0p .* Z0pp) ./ s_theta .^ 2 .* d_s_thetatheta_dc(:,mode);
    end
end
function[d_k_dc] = compute_d_k_dc(s_theta, s_thetatheta,  ...
                                  d_s_theta_dc, d_s_thetatheta_dc,...
                                  r_theta_cyl, r_thetatheta_cyl, ...
                                  d_r_theta_cyl_dc, d_r_thetatheta_cyl_dc,...
                                  phi, R0c)

    d_k_dc = zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d_k_dc(:,mode) = (d_r_thetatheta_cyl_dc(:,mode) .* s_theta .^ 2 + (-2 .* d_s_theta_dc(:,mode) .* r_thetatheta_cyl - r_theta_cyl .* d_s_thetatheta_dc(:,mode) - s_thetatheta .* d_r_theta_cyl_dc(:,mode)) .* s_theta + 3 .* r_theta_cyl .* s_thetatheta .* d_s_theta_dc(:,mode)) ./ s_theta .^ 4;
    end
end
function [d_k_dtheta_dc] = d_kappa_dtheta_dc(s_theta, s_thetatheta, s_thetathetatheta, ...
                                             d_s_theta_dc, d_s_thetatheta_dc,d_s_thetathetatheta_dc,...
                                             r_theta_cyl, r_thetatheta_cyl, r_theta_cyl_tc, r_thetatheta_cyl_tc,...
                                             d_r_theta_cyl_dc, d_r_thetatheta_cyl_dc,d_r_theta_cyl_tc_dc, d_r_thetatheta_cyl_tc_dc,...
                                             phi, R0c)
    d_k_dtheta_dc = zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d_k_dtheta_dc(:,mode) =  (d_r_thetatheta_cyl_tc_dc(:,mode) .* s_theta .^ 3 + ((-d_r_theta_cyl_tc_dc(:,mode) - 2 .* d_r_thetatheta_cyl_dc(:,mode)) .* s_thetatheta - 2 .* r_thetatheta_cyl_tc .* d_s_theta_dc(:,mode) - d_s_thetathetatheta_dc(:,mode) .* r_theta_cyl + (-r_theta_cyl_tc - 2 .* r_thetatheta_cyl) .* d_s_thetatheta_dc(:,mode) - d_r_theta_cyl_dc(:,mode) .* s_thetathetatheta) .* s_theta .^ 2 + (3 .* d_r_theta_cyl_dc(:,mode) .* s_thetatheta .^ 2 + ((3 .* r_theta_cyl_tc + 6 .* r_thetatheta_cyl) .* d_s_theta_dc(:,mode) + 6 .* r_theta_cyl .* d_s_thetatheta_dc(:,mode)) .* s_thetatheta + 3 .* r_theta_cyl .* d_s_theta_dc(:,mode) .* s_thetathetatheta) .* s_theta - 12 .* d_s_theta_dc(:,mode) .* r_theta_cyl .* s_thetatheta .^ 2) ./ s_theta .^ 5;
    end
end
function[d_curvature_dt_dc] = compute_dk_dt_dc(k1, k2, k3,...
                                               d_k1_dtheta, d_k2_dtheta, d_k3_dtheta,...
                                               d_k1_dc, d_k2_dc, d_k3_dc,...
                                               d_k1_dtheta_dc, d_k2_dtheta_dc, d_k3_dtheta_dc,...
                                               curvature,d_curvature_dc,...
                                               phi, R0c)
    d_curvature_dt_dc = zeros(size(phi,1), numel(R0c));

    for mode = 1:numel(R0c)
        d_curvature_dt_dc(:,mode) = (d_k1_dc(:,mode) .* d_k1_dtheta + d_k2_dc(:,mode) .* d_k2_dtheta + d_k3_dc(:,mode) .* d_k3_dtheta + k1 .* d_k1_dtheta_dc(:,mode) + k2 .* d_k2_dtheta_dc(:,mode) + k3 .* d_k3_dtheta_dc(:,mode)) ./ curvature ...
                                  - (k1 .* d_k1_dtheta + k2 .* d_k2_dtheta + k3 .* d_k3_dtheta) ./ curvature .^ 2 .* d_curvature_dc(:,mode);
    end

end
function [d_torsion_dc] = compute_torsion_dc( torsion_num, torsion_denom,...
                                              xp, xpp, xppp, d_xp_dc, d_xpp_dc, d_xppp_dc, ...
                                              yp, ypp, yppp, d_yp_dc, d_ypp_dc, d_yppp_dc, ...
                                              zp, zpp, zppp, d_zp_dc, d_zpp_dc, d_zppp_dc,...
                                              phi, R0c)
    d_torsion_dc = zeros(size(torsion_num,1), numel(R0c));

    for mode = 1:numel(R0c)
        d_torsion_num_dc = d_xppp_dc(:,mode) .* (yp .* zpp - ypp .* zp) + xppp .* (d_yp_dc(:,mode) .* zpp - ypp .* d_zp_dc(:,mode) + yp .* d_zpp_dc(:,mode) - d_ypp_dc(:,mode) .* zp) ...
                         + d_yppp_dc(:,mode) .* (-xp .* zpp + xpp .* zp)+ yppp .* (-d_xp_dc(:,mode) .* zpp + xpp .* d_zp_dc(:,mode) - xp .* d_zpp_dc(:,mode) + d_xpp_dc(:,mode) .* zp) ...
                         + d_zppp_dc(:,mode) .* (xp .* ypp - xpp .* yp)  + zppp .* (d_xp_dc(:,mode) .* ypp - xpp .* d_yp_dc(:,mode) + xp .* d_ypp_dc(:,mode) - d_xpp_dc(:,mode) .* yp);
        d_torsion_denom_dc = 2 * (yp .* zpp - ypp .* zp) .* (d_yp_dc(:,mode) .* zpp - ypp .* d_zp_dc(:,mode) + yp .* d_zpp_dc(:,mode) - d_ypp_dc(:,mode) .* zp) ...
                           + 2 .* (-xp .* zpp + xpp .* zp) .* (-d_xp_dc(:,mode) .* zpp + xpp .* d_zp_dc(:,mode) - xp .* d_zpp_dc(:,mode) + d_xpp_dc(:,mode) .* zp) ...
                           + 2 * (xp .* ypp - xpp .* yp) .* (d_xp_dc(:,mode) .* ypp - xpp .* d_yp_dc(:,mode) + xp .* d_ypp_dc(:,mode) - d_xpp_dc(:,mode) .* yp);


        d_torsion_dc(:,mode) = (d_torsion_num_dc .* torsion_denom - d_torsion_denom_dc .* torsion_num)./torsion_denom.^2;

        for i = 1:size(phi,1)
            if( abs(d_torsion_num_dc(i)) == 0. && abs(d_torsion_denom_dc(i)) == 0. ) 
                d_torsion_dc(i,mode) = 0.;
            end
        end

    end

end


function [d2_X1c_dc] = compute_d2_X1c_dc(abs_G_over_B0, d_abs_G_over_B0_dc,...
                                         s_theta, d_s_theta_dc,...
                                         curvature, d_curvature_dc, d_curvature_dtheta, d_curvature_dtheta_dc,...
                                         phi, R0c, eta_bar)
    d2_X1c_dc =  zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d2_X1c_dc(:,mode) = -eta_bar .* (((abs_G_over_B0 .* d_curvature_dtheta_dc(:,mode) + d_abs_G_over_B0_dc(:,mode) .* d_curvature_dtheta) .* s_theta - d_s_theta_dc(:,mode) .* abs_G_over_B0 .* d_curvature_dtheta) .* curvature - 2 .* s_theta .* abs_G_over_B0 .* d_curvature_dc(:,mode) .* d_curvature_dtheta) ./ s_theta .^ 2 ./ curvature .^ 3;
    end
end

function [d2_Y1s_dc] = compute_d2_Y1s_dc(abs_G_over_B0, d_abs_G_over_B0_dc,...
                                         s_theta, d_s_theta_dc,...
                                         curvature, d_curvature_dtheta, d_curvature_dtheta_dc,...
                                         phi, R0c, eta_bar)
    d2_Y1s_dc =  zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d2_Y1s_dc(:,mode) = (abs_G_over_B0 .* d_curvature_dtheta_dc(:,mode) .* s_theta + d_abs_G_over_B0_dc(:,mode) .* d_curvature_dtheta .* s_theta ...
                           - abs_G_over_B0 .* d_curvature_dtheta .* d_s_theta_dc(:,mode)) ./ s_theta .^ 2 ./ eta_bar;
    end
end

function [d2_Y1c_dc] = compute_d2_Y1c_dc(abs_G_over_B0, d_abs_G_over_B0_dc,...
                                         s_theta, d_s_theta_dc,...
                                         curvature, d_curvature_dc, d_curvature_dtheta, d_curvature_dtheta_dc,...
                                         phi, R0c, eta_bar, D, sigma)
    d2_Y1c_dc =  zeros(size(phi,1), numel(R0c));
    for mode = 1:numel(R0c)
        d2_Y1c_dc(:,mode) = (s_theta .* (D*sigma .* curvature + d_curvature_dtheta .* sigma) .* d_abs_G_over_B0_dc(:,mode) ...
                          + ((-D*sigma .* curvature - d_curvature_dtheta .* sigma) .* d_s_theta_dc(:,mode) ...
                          + s_theta .* (D*sigma .* d_curvature_dc(:,mode) + d_curvature_dtheta_dc(:,mode) .* sigma)) .* abs_G_over_B0) ./ s_theta .^ 2 ./ eta_bar;
    end
end
function [d_B1x_dc, d_B1y_dc, d_B1z_dc,...
          d_B2x_dc, d_B2y_dc, d_B2z_dc,...
          d_B3x_dc, d_B3y_dc, d_B3z_dc] = compute_gradB_dc(abs_G_over_B0,d_abs_G0_over_B0_dc,...
                                                           torsion, curvature, ...
                                                           d_torsion_dc, d_curvature_dc,...
                                                           t,n,b,...
                                                           d_t1_dc, d_n1_dc,d_b1_dc,...
                                                           d_t2_dc, d_n2_dc,d_b2_dc,...
                                                           d_t3_dc, d_n3_dc,d_b3_dc,...
                                                           X1c, Y1s, Y1c,d_X1c_dc, d_Y1s_dc, d_Y1c_dc,...
                                                           dX1c, dY1s, dY1c, d2_X1c_dc,d2_Y1s_dc,d2_Y1c_dc,...
                                                           s_theta, d_stheta_dc, B0, iota, phi, R0c)
    svp = 1;
    sg = 1;

    d_B1x_dc = zeros(size(phi,1),numel(R0c));
    d_B1y_dc = zeros(size(phi,1),numel(R0c));
    d_B1z_dc = zeros(size(phi,1),numel(R0c));

    tj = t(:,1);
    nj = n(:,1);
    bj = b(:,1);

    d_tj_dc = d_t1_dc;
    d_nj_dc = d_n1_dc;
    d_bj_dc = d_b1_dc;

    for mode = 1:numel(R0c)
        d_n_dc = [d_n1_dc(:,mode) d_n2_dc(:,mode) d_n3_dc(:,mode)];
        d_t_dc = [d_t1_dc(:,mode) d_t2_dc(:,mode) d_t3_dc(:,mode)];
        d_b_dc = [d_b1_dc(:,mode) d_b2_dc(:,mode) d_b3_dc(:,mode)];

        out =   -(abs_G_over_B0 ^ 2 * ((-bj .* d_torsion_dc(:,mode) - d_bj_dc(:,mode) .* torsion - d_tj_dc(:,mode) .* curvature - tj .* d_curvature_dc(:,mode)) .* n - d_n_dc .* torsion .* bj + (b .* d_torsion_dc(:,mode) + d_b_dc .* torsion) .* nj + b .* d_nj_dc(:,mode) .* torsion - curvature .* tj .* d_n_dc) * svp^ 2 + (((((-2 * Y1c .* d_Y1c_dc(:,mode) - 2 * Y1s .* d_Y1s_dc(:,mode)) * iota + d2_Y1s_dc(:,mode) .* Y1c - d2_Y1c_dc(:,mode) .* Y1s - dY1c .* d_Y1s_dc(:,mode) + dY1s .* d_Y1c_dc(:,mode)) .* bj + ((-Y1c .* d_nj_dc(:,mode) - nj .* d_Y1c_dc(:,mode)) .* X1c - Y1c.^ 2 .* d_bj_dc(:,mode) - d_X1c_dc(:,mode) .* nj.* Y1c - Y1s.^ 2 .* d_bj_dc(:,mode)) * iota + Y1c.* dY1s .* d_bj_dc(:,mode) + (-Y1s .* d2_X1c_dc(:,mode) - dX1c .* d_Y1s_dc(:,mode)) .* nj - Y1s .* (dX1c .* d_nj_dc(:,mode) + dY1c .* d_bj_dc(:,mode))) .* n ...
                + (((Y1c .* d_b_dc + b .* d_Y1c_dc(:,mode)) .* X1c + b .* Y1c .* d_X1c_dc(:,mode) - d_n_dc .* Y1c .^ 2 - d_n_dc .* Y1s .^ 2) .* iota + (-b .* d2_Y1s_dc(:,mode) - dY1s .* d_b_dc) .* X1c - b .* dY1s .* d_X1c_dc(:,mode) + dY1s .* d_n_dc .* Y1c - d_n_dc .* dY1c .* Y1s) .* bj + X1c .* ((b .* d_nj_dc(:,mode) + d_b_dc .* nj) .* X1c + (b .* d_bj_dc(:,mode) - d_n_dc .* nj) .* Y1c + 2 .* b .* d_X1c_dc(:,mode) .* nj) .* iota - b .* dY1s .* d_bj_dc(:,mode) .* X1c - Y1s .* dX1c .* d_n_dc .* nj) .* abs_G_over_B0 - ((((-Y1c .^ 2 - Y1s .^ 2) .* iota + dY1s .* Y1c - dY1c .* Y1s) .* bj - (X1c .* Y1c .* iota + Y1s .* dX1c) .* nj) .* n + X1c .* ((Y1c .* iota - dY1s) .* bj + iota .* nj .* X1c) .* b) .* d_abs_G0_over_B0_dc(:,mode)) .* svp - abs_G_over_B0 .^ 2 .* ((d_t_dc .* curvature + t .* d_curvature_dc(:,mode)) .* nj + d_nj_dc(:,mode) .* curvature .* t) .* sg) .* B0 ./ abs_G_over_B0 .^ 2;

        d_B1x_dc(:,mode) = out(:,1);
        d_B1y_dc(:,mode) = out(:,2);
        d_B1z_dc(:,mode) = out(:,3);
    end








    d_B2x_dc = zeros(size(phi,1),numel(R0c));
    d_B2y_dc = zeros(size(phi,1),numel(R0c));
    d_B2z_dc = zeros(size(phi,1),numel(R0c));

    tj = t(:,2);
    nj = n(:,2);
    bj = b(:,2);

    d_tj_dc = d_t2_dc;
    d_nj_dc = d_n2_dc;
    d_bj_dc = d_b2_dc;


    for mode = 1:numel(R0c)
        d_n_dc = [d_n1_dc(:,mode) d_n2_dc(:,mode) d_n3_dc(:,mode)];
        d_t_dc = [d_t1_dc(:,mode) d_t2_dc(:,mode) d_t3_dc(:,mode)];
        d_b_dc = [d_b1_dc(:,mode) d_b2_dc(:,mode) d_b3_dc(:,mode)];
        out = -(abs_G_over_B0 .^ 2 .* ((-bj .* d_torsion_dc(:,mode) - d_bj_dc(:,mode) .* torsion - d_tj_dc(:,mode) .* curvature - tj .* d_curvature_dc(:,mode)) .* n - d_n_dc .* torsion .* bj + (b .* d_torsion_dc(:,mode) + d_b_dc .* torsion) .* nj + b .* d_nj_dc(:,mode) .* torsion - curvature .* tj .* d_n_dc) .* svp .^ 2 + (((((-2 .* Y1c .* d_Y1c_dc(:,mode) - 2 .* Y1s .* d_Y1s_dc(:,mode)) .* iota + d2_Y1s_dc(:,mode) .* Y1c - d2_Y1c_dc(:,mode) .* Y1s - dY1c .* d_Y1s_dc(:,mode) + dY1s .* d_Y1c_dc(:,mode)) .* bj + ((-Y1c .* d_nj_dc(:,mode) - nj .* d_Y1c_dc(:,mode)) .* X1c - Y1c .^ 2 .* d_bj_dc(:,mode) - d_X1c_dc(:,mode) .* nj .* Y1c - Y1s .^ 2 .* d_bj_dc(:,mode)) .* iota + Y1c .* dY1s .* d_bj_dc(:,mode) + (-Y1s .* d2_X1c_dc(:,mode) - dX1c .* d_Y1s_dc(:,mode)) .* nj - Y1s .* (dX1c .* d_nj_dc(:,mode) + dY1c .* d_bj_dc(:,mode))) .* n + (((Y1c .* d_b_dc + b .* d_Y1c_dc(:,mode)) .* X1c + b .* Y1c .* d_X1c_dc(:,mode) - d_n_dc .* Y1c .^ 2 - d_n_dc .* Y1s .^ 2) .* iota + (-b .* d2_Y1s_dc(:,mode) - dY1s .* d_b_dc) .* X1c - b .* dY1s .* d_X1c_dc(:,mode) + dY1s .* d_n_dc .* Y1c - d_n_dc .* dY1c .* Y1s) .* bj + X1c .* ((b .* d_nj_dc(:,mode) + d_b_dc .* nj) .* X1c + (b .* d_bj_dc(:,mode) - d_n_dc .* nj) .* Y1c + 2 .* b .* d_X1c_dc(:,mode) .* nj) .* iota - b .* dY1s .* d_bj_dc(:,mode) .* X1c - Y1s .* dX1c .* d_n_dc .* nj) .* abs_G_over_B0 - ((((-Y1c .^ 2 - Y1s .^ 2) .* iota + dY1s .* Y1c - dY1c .* Y1s) .* bj - (X1c .* Y1c .* iota + Y1s .* dX1c) .* nj) .* n + X1c .* ((Y1c .* iota - dY1s) .* bj + iota .* nj .* X1c) .* b) .* d_abs_G0_over_B0_dc(:,mode)) .* svp - abs_G_over_B0 .^ 2 .* ((d_t_dc .* curvature + t .* d_curvature_dc(:,mode)) .* nj + d_nj_dc(:,mode) .* curvature .* t) .* sg) .* B0 ./ abs_G_over_B0 .^ 2;
        d_B2x_dc(:,mode) = out(:,1);
        d_B2y_dc(:,mode) = out(:,2);
        d_B2z_dc(:,mode) = out(:,3);
    end








    d_B3x_dc = zeros(size(phi,1),numel(R0c));
    d_B3y_dc = zeros(size(phi,1),numel(R0c));
    d_B3z_dc = zeros(size(phi,1),numel(R0c));

    tj = t(:,3);
    nj = n(:,3);
    bj = b(:,3);

    d_tj_dc = d_t3_dc;
    d_nj_dc = d_n3_dc;
    d_bj_dc = d_b3_dc;


    for mode = 1:numel(R0c)
        d_n_dc = [d_n1_dc(:,mode) d_n2_dc(:,mode) d_n3_dc(:,mode)];
        d_t_dc = [d_t1_dc(:,mode) d_t2_dc(:,mode) d_t3_dc(:,mode)];
        d_b_dc = [d_b1_dc(:,mode) d_b2_dc(:,mode) d_b3_dc(:,mode)];
        out = -(abs_G_over_B0 .^ 2 .* ((-bj .* d_torsion_dc(:,mode) - d_bj_dc(:,mode) .* torsion - d_tj_dc(:,mode) .* curvature - tj .* d_curvature_dc(:,mode)) .* n - d_n_dc .* torsion .* bj + (b .* d_torsion_dc(:,mode) + d_b_dc .* torsion) .* nj + b .* d_nj_dc(:,mode) .* torsion - curvature .* tj .* d_n_dc) .* svp .^ 2 + (((((-2 .* Y1c .* d_Y1c_dc(:,mode) - 2 .* Y1s .* d_Y1s_dc(:,mode)) .* iota + d2_Y1s_dc(:,mode) .* Y1c - d2_Y1c_dc(:,mode) .* Y1s - dY1c .* d_Y1s_dc(:,mode) + dY1s .* d_Y1c_dc(:,mode)) .* bj + ((-Y1c .* d_nj_dc(:,mode) - nj .* d_Y1c_dc(:,mode)) .* X1c - Y1c .^ 2 .* d_bj_dc(:,mode) - d_X1c_dc(:,mode) .* nj .* Y1c - Y1s .^ 2 .* d_bj_dc(:,mode)) .* iota + Y1c .* dY1s .* d_bj_dc(:,mode) + (-Y1s .* d2_X1c_dc(:,mode) - dX1c .* d_Y1s_dc(:,mode)) .* nj - Y1s .* (dX1c .* d_nj_dc(:,mode) + dY1c .* d_bj_dc(:,mode))) .* n + (((Y1c .* d_b_dc + b .* d_Y1c_dc(:,mode)) .* X1c + b .* Y1c .* d_X1c_dc(:,mode) - d_n_dc .* Y1c .^ 2 - d_n_dc .* Y1s .^ 2) .* iota + (-b .* d2_Y1s_dc(:,mode) - dY1s .* d_b_dc) .* X1c - b .* dY1s .* d_X1c_dc(:,mode) + dY1s .* d_n_dc .* Y1c - d_n_dc .* dY1c .* Y1s) .* bj + X1c .* ((b .* d_nj_dc(:,mode) + d_b_dc .* nj) .* X1c + (b .* d_bj_dc(:,mode) - d_n_dc .* nj) .* Y1c + 2 .* b .* d_X1c_dc(:,mode) .* nj) .* iota - b .* dY1s .* d_bj_dc(:,mode) .* X1c - Y1s .* dX1c .* d_n_dc .* nj) .* abs_G_over_B0 - ((((-Y1c .^ 2 - Y1s .^ 2) .* iota + dY1s .* Y1c - dY1c .* Y1s) .* bj - (X1c .* Y1c .* iota + Y1s .* dX1c) .* nj) .* n + X1c .* ((Y1c .* iota - dY1s) .* bj + iota .* nj .* X1c) .* b) .* d_abs_G0_over_B0_dc(:,mode)) .* svp - abs_G_over_B0 .^ 2 .* ((d_t_dc .* curvature + t .* d_curvature_dc(:,mode)) .* nj + d_nj_dc(:,mode) .* curvature .* t) .* sg) .* B0 ./ abs_G_over_B0 .^ 2;
        d_B3x_dc(:,mode) = out(:,1);
        d_B3y_dc(:,mode) = out(:,2);
        d_B3z_dc(:,mode) = out(:,3);
    end


end



function [d_gradB1_detabar, d_gradB2_detabar, d_gradB3_detabar] = compute_gradB_detabar(B0, abs_G_over_B0,...
                                                                                        iota, ...
                                                                                        t,n,b,...
                                                                                        X1c, Y1s, Y1c, dX1c, dY1s, dY1c,...
                                                                                        d_X1c_detabar, d_Y1s_detabar, d_Y1c_detabar, ...
                                                                                        d2_X1c_detabar, d2_Y1s_detabar, d2_Y1c_detabar)
        svp = 1;
        sg = 1;

        tj = t(:,1);
        nj = n(:,1);
        bj = b(:,1);
        d_gradB1_detabar = svp .* ((((d2_Y1c_detabar .* Y1s + dY1c .* d_Y1s_detabar - d2_Y1s_detabar .* Y1c - dY1s .* d_Y1c_detabar + iota .* (2 .* Y1c .* d_Y1c_detabar + 2 .* Y1s .* d_Y1s_detabar)) .* n + ((-d_Y1c_detabar .* X1c - d_X1c_detabar .* Y1c) .* iota + d_X1c_detabar .* dY1s + X1c .* d2_Y1s_detabar) .* b) .* bj) - 0.2e1 .* nj .* (((-(d_X1c_detabar .* Y1c) ./ 0.2e1 - (d_Y1c_detabar .* X1c) ./ 0.2e1) .* iota - (dX1c .* d_Y1s_detabar) ./ 0.2e1 - (d2_X1c_detabar .* Y1s) ./ 0.2e1) .* n + (b .* iota .* d_X1c_detabar .* X1c))) .* B0 ./ abs_G_over_B0;
        
        tj = t(:,2);
        nj = n(:,2);
        bj = b(:,2);
        d_gradB2_detabar = svp .* ((((d2_Y1c_detabar .* Y1s + dY1c .* d_Y1s_detabar - d2_Y1s_detabar .* Y1c - dY1s .* d_Y1c_detabar + iota .* (2 .* Y1c .* d_Y1c_detabar + 2 .* Y1s .* d_Y1s_detabar)) .* n + ((-d_Y1c_detabar .* X1c - d_X1c_detabar .* Y1c) .* iota + d_X1c_detabar .* dY1s + X1c .* d2_Y1s_detabar) .* b) .* bj) - 0.2e1 .* nj .* (((-(d_X1c_detabar .* Y1c) ./ 0.2e1 - (d_Y1c_detabar .* X1c) ./ 0.2e1) .* iota - (dX1c .* d_Y1s_detabar) ./ 0.2e1 - (d2_X1c_detabar .* Y1s) ./ 0.2e1) .* n + (b .* iota .* d_X1c_detabar .* X1c))) .* B0 ./ abs_G_over_B0;
        
        tj = t(:,3);
        nj = n(:,3);
        bj = b(:,3);
        d_gradB3_detabar = svp .* ((((d2_Y1c_detabar .* Y1s + dY1c .* d_Y1s_detabar - d2_Y1s_detabar .* Y1c - dY1s .* d_Y1c_detabar + iota .* (2 .* Y1c .* d_Y1c_detabar + 2 .* Y1s .* d_Y1s_detabar)) .* n + ((-d_Y1c_detabar .* X1c - d_X1c_detabar .* Y1c) .* iota + d_X1c_detabar .* dY1s + X1c .* d2_Y1s_detabar) .* b) .* bj) - 0.2e1 .* nj .* (((-(d_X1c_detabar .* Y1c) ./ 0.2e1 - (d_Y1c_detabar .* X1c) ./ 0.2e1) .* iota - (dX1c .* d_Y1s_detabar) ./ 0.2e1 - (d2_X1c_detabar .* Y1s) ./ 0.2e1) .* n + (b .* iota .* d_X1c_detabar .* X1c))) .* B0 ./ abs_G_over_B0;
        
end









function [d_B1x_ds, d_B1y_ds, d_B1z_ds,...
         d_B2x_ds, d_B2y_ds, d_B2z_ds,...
         d_B3x_ds, d_B3y_ds, d_B3z_ds] = compute_gradB_ds(B0, abs_G_over_B0,...
                                                          phi, sigma,iota,...
                                                          n,b,...
                                                          X1c, Y1s, Y1c,d_Y1c_ds,dY1s, d2_Y1c_ds)
    svp = 1;

    d_B1x_ds = zeros(size(phi,1),size(sigma,1));
    d_B1y_ds = zeros(size(phi,1),size(sigma,1));
    d_B1z_ds = zeros(size(phi,1),size(sigma,1));

    nj = n(:,1);
    bj = b(:,1);

    for pt = 1:numel(phi)
        d_B1x_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,1) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,1));
        d_B1y_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,2) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,2));
        d_B1z_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,3) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,3));
    end

    d_B2x_ds = zeros(size(phi,1),size(sigma,1));
    d_B2y_ds = zeros(size(phi,1),size(sigma,1));
    d_B2z_ds = zeros(size(phi,1),size(sigma,1));

    nj = n(:,2);
    bj = b(:,2);

    for pt = 1:numel(phi)
        d_B2x_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,1) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,1));
        d_B2y_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,2) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,2));
        d_B2z_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,3) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,3));
    end


    d_B3x_ds = zeros(size(phi,1),size(sigma,1));
    d_B3y_ds = zeros(size(phi,1),size(sigma,1));
    d_B3z_ds = zeros(size(phi,1),size(sigma,1));

    nj = n(:,3);
    bj = b(:,3);

    for pt = 1:numel(phi)
        d_B3x_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,1) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,1));
        d_B3y_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,2) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,2));
        d_B3z_ds(pt,:)  = svp * B0 / abs_G_over_B0 * ((iota * X1c(pt) * d_Y1c_ds(pt,:) * nj(pt) + (2 * iota * Y1c(pt) .* d_Y1c_ds(pt,:) + d2_Y1c_ds(pt,:) * Y1s(pt) - dY1s(pt) * d_Y1c_ds(pt,:)) * bj(pt)) .* n(pt,3) - iota .* X1c(pt) .* d_Y1c_ds(pt,:) * bj(pt) * b(pt,3));
    end

end

function [d_gradB1_diota,d_gradB2_diota,d_gradB3_diota] = compute_gradB_diota(B0, abs_G_over_B0,...
                                                                  n,b,X1c, Y1s, Y1c)
    svp = 1;


    nj = n(:,1);
    bj = b(:,1);

    d_gradB1_diota =   svp * B0 / abs_G_over_B0 * ((X1c .* Y1c .* nj + (Y1c .^ 2 + Y1s .^ 2) .* bj) .* n + (-X1c .^ 2 .* nj - X1c .* Y1c .* bj) .* b);

    nj = n(:,2);
    bj = b(:,2);
    d_gradB2_diota =   svp * B0 / abs_G_over_B0 * ((X1c .* Y1c .* nj + (Y1c .^ 2 + Y1s .^ 2) .* bj) .* n + (-X1c .^ 2 .* nj - X1c .* Y1c .* bj) .* b);

    nj = n(:,3);
    bj = b(:,3);
    d_gradB3_diota =   svp * B0 / abs_G_over_B0 * ((X1c .* Y1c .* nj + (Y1c .^ 2 + Y1s .^ 2) .* bj) .* n + (-X1c .^ 2 .* nj - X1c .* Y1c .* bj) .* b);

end
        

