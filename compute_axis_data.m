function [BR,Bphi,BZ,...
          d_BR_d_R,d_Bphi_d_R,d_BZ_d_R,...
          d_BR_d_phi,d_Bphi_d_phi,d_BZ_d_phi,...
          d_BR_d_Z,d_Bphi_d_Z,d_BZ_d_Z,...
          dB_R_dc, dB_phi_dc, dB_Z_dc, ...
          d2_R_BR_Bphi_dR_dc, d2_R_BR_Bphi_dZ_dc, d2_R_BZ_Bphi_dR_dc, d2_R_BZ_Bphi_dZ_dc,...
          d2_R_BR_Bphi_d2R, d2_R_BR_Bphi_dZ_dR, d2_R_BZ_Bphi_d2R, d2_R_BZ_Bphi_dZ_dR,...
          d2_R_BR_Bphi_d2Z, d2_R_BZ_Bphi_d2Z,...
          dBR_dI, dBphi_dI, dBZ_dI,...
          d2_R_BR_Bphi_dR_dI, d2_R_BR_Bphi_dZ_dI, ...
          d2_R_BZ_Bphi_dR_dI, d2_R_BZ_Bphi_dZ_dI] = compute_axis_data(R,phi,Z, coilData, fourierData)

    x = R * cos(phi);
    y = R * sin(phi);
    z = Z;

    x_R = [x/R, y/R, 0];
    x_Theta= [-y, x, 0];
    x_Z = [0, 0, 1];



    %% Declaration of variables
    gamma = 4 * pi * 10^(-7) ;
    dB = zeros(size(coilData.coil_tangents{1},1),3);
    dB_dR = zeros(size(coilData.coil_tangents{1},1),3);
    dB_dTheta = zeros(size(coilData.coil_tangents{1},1),3);
    dB_dZ = zeros(size(coilData.coil_tangents{1},1),3);
    dB_xyz_dc = zeros(numel(coilData.coil_coeffs),3);

    d2_B_d2R = zeros(size(coilData.coil_tangents{1},1),3);
    d2_B_dR_dZ = zeros(size(coilData.coil_tangents{1},1),3);
    d2_B_d2Z = zeros(size(coilData.coil_tangents{1},1),3);
    
    d2_Bxyz_dR_dc = zeros(numel(coilData.coil_coeffs),3);
    d2_Bxyz_dZ_dc = zeros(numel(coilData.coil_coeffs),3);
    
    dB_dI = zeros(size(coilData.I,1),3);
    d2B_dR_dI = zeros(size(coilData.I,1),3);
    d2B_dZ_dI = zeros(size(coilData.I,1),3);
    
    point = [x y z];
    for i = 1: coilData.C
        %% Numerical integration of Biot-Savart law
        diff = point-coilData.coil_field{i};
        d_diff_dR = repmat(x_R, size(coilData.coil_tangents{i},1), 1);
        d_diff_dZ = repmat(x_Z, size(coilData.coil_tangents{i},1), 1);
        
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
        dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
        dist_Z     = diff(:,3)./dist;

        dist3_R        = 3 * dist.^2 .* dist_R;
        dist3_Theta    = 3 * dist.^2 .* dist_Theta;
        dist3_Z        = 3 * dist.^2 .* dist_Z;
        distm4_R       =-4 * dist.^(-5) .* dist_R;
        distm4_Z       =-4 * dist.^(-5) .* dist_Z;
        d2_dist_dR2    = (dist - dist_R.* ((x/R) * diff(:,1) + (y/R) * diff(:,2))    )./(dist.^2);
        d2_dist_dR_dZ  = -diff(:,3)./dist.^2 .* dist_R;
        d2_dist_d2Z    = ( dist - dist_Z.*diff(:,3) )./ dist.^2;

        T_x_diff = cross(coilData.coil_tangents{i}, diff);
        T_x_diff_R = cross(coilData.coil_tangents{i}, repmat(x_R, size(coilData.coil_tangents{i},1), 1));
        T_x_diff_Theta = cross(coilData.coil_tangents{i}, repmat(x_Theta, size(coilData.coil_tangents{i},1), 1));
        T_x_diff_Z = cross(coilData.coil_tangents{i}, repmat(x_Z, size(coilData.coil_tangents{i},1), 1));

        
        term1RR =  - (3./dist.^4) .* dist_R .* T_x_diff_R;
        term2RR =  distm4_R .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR2 .* T_x_diff + T_x_diff_R .* dist_R);
        d2_B_d2R = d2_B_d2R + (gamma * coilData.I(i) ) / (4. * pi) * (term1RR-3.*term2RR);

        term1RZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_R;
        term2RZ =  distm4_Z .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR_dZ .* T_x_diff + T_x_diff_Z .* dist_R);
        d2_B_dR_dZ = d2_B_dR_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (term1RZ-3.*term2RZ);

        term1ZZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_Z;
        term2ZZ =  distm4_Z .* dist_Z .* T_x_diff + (1./dist.^4).*(d2_dist_d2Z .* T_x_diff + T_x_diff_Z .* dist_Z);
        d2_B_d2Z = d2_B_d2Z + (gamma * coilData.I(i) ) / (4. * pi) * (term1ZZ-3.*term2ZZ);
        
        dB = dB + (gamma * coilData.I(i) ) / (4. * pi) * T_x_diff./dist.^3;
        dB_dR = dB_dR + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
        dB_dTheta = dB_dTheta + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
        dB_dZ = dB_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;
        
        dB_dI(i,:) = dB_dI(i,:) + sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
        d2B_dR_dI(i,:) = d2B_dR_dI(i,:) + sum(1. / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.,1);
        d2B_dZ_dI(i,:) = d2B_dZ_dI(i,:) + sum(1. / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.,1);
        
        for coord = 1:3
            for sc = 1:2 % sine then cos
                for k = 1:coilData.Nt
                    idx = (i - 1) * coilData.Nt * 3 * 2 + (coord - 1) * coilData.Nt * 2 + (sc-1)*coilData.Nt + k;

                    dT_dc = zeros(coilData.M,3);
                    dxyz_dc = zeros(coilData.M,3);

                    if sc == 1
                        dxyz_dc(:,coord) = fourierData.SIN(:,k);
                        dT_dc(:,coord) = fourierData.dSIN(:,k);
                    else
                        dxyz_dc(:,coord) = fourierData.COS(:,k);
                        dT_dc(:,coord) = fourierData.dCOS(:,k);
                    end
                    d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                    d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;
                    d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                    
                    d2_dist_dR_dc = dot(repmat(x_R, size(coilData.coil_tangents{i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);
                    d2_dist_dZ_dc = dot(repmat(x_Z, size(coilData.coil_tangents{i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);

                    d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{i}, -dxyz_dc);
                    dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);

                    term1 =  cross(dT_dc, d_diff_dR) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{i}, d_diff_dR)./dist.^6;
                    term2 = d_distm4_dc .* dist_R .* T_x_diff + (1./dist).^4 .* (d2_dist_dR_dc .* T_x_diff + d_T_x_diff_dc .* dist_R );

                    term1z =  cross(dT_dc, d_diff_dZ) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{i}, d_diff_dZ)./dist.^6;
                    term2z = d_distm4_dc .* dist_Z .* T_x_diff + (1./dist).^4 .* (d2_dist_dZ_dc .* T_x_diff + d_T_x_diff_dc .* dist_Z );

                    d2_Bxyz_dR_dc(idx,:) = d2_Bxyz_dR_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(term1 - 3 * term2,1);
                    d2_Bxyz_dZ_dc(idx,:) = d2_Bxyz_dZ_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(term1z - 3 * term2z,1);
                end
            end
        end


        if coilData.ss == 1
            diff = point-coilData.coil_field{coilData.ts*coilData.C + i};
            d_diff_dR = repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1);
            d_diff_dZ = repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1);
            
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
            dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
            dist_Z     = diff(:,3)./dist;

            dist3_R        = 3 * dist.^2 .* dist_R;
            dist3_Theta    = 3 * dist.^2 .* dist_Theta;
            dist3_Z        = 3 * dist.^2 .* dist_Z; 
            distm4_R       =-4 * dist.^(-5) .* dist_R;
            distm4_Z       =-4 * dist.^(-5) .* dist_Z;
            d2_dist_dR2    = (dist - dist_R.* ((x/R) * diff(:,1) + (y/R) * diff(:,2))    )./(dist.^2);
            d2_dist_dR_dZ  = -diff(:,3)./dist.^2 .* dist_R;
            d2_dist_d2Z    = ( dist - dist_Z.*diff(:,3) )./ dist.^2;
            
            T_x_diff = cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, diff);
            T_x_diff_R = cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1));
            T_x_diff_Theta = cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1));
            T_x_diff_Z = cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1));

            
            term1RR =  - (3./dist.^4) .* dist_R .* T_x_diff_R;
            term2RR =  distm4_R .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR2 .* T_x_diff + T_x_diff_R .* dist_R);
            d2_B_d2R = d2_B_d2R + (gamma * -coilData.I(i) ) / (4. * pi) * (term1RR-3.*term2RR);

            term1RZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_R;
            term2RZ =  distm4_Z .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR_dZ .* T_x_diff + T_x_diff_Z .* dist_R);
            d2_B_dR_dZ = d2_B_dR_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (term1RZ-3.*term2RZ);
            
            term1ZZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_Z;
            term2ZZ =  distm4_Z .* dist_Z .* T_x_diff + (1./dist.^4).*(d2_dist_d2Z .* T_x_diff + T_x_diff_Z .* dist_Z);
            d2_B_d2Z = d2_B_d2Z + (gamma * -coilData.I(i) ) / (4. * pi) * (term1ZZ-3.*term2ZZ);
            
            dB = dB + (gamma * -coilData.I(i) ) / (4. * pi) * T_x_diff./dist.^3;
            dB_dR = dB_dR + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
            dB_dTheta = dB_dTheta + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
            dB_dZ = dB_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;

            dB_dI(i,:) = dB_dI(i,:) - sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
            d2B_dR_dI(i,:) = d2B_dR_dI(i,:) - sum(1. / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.,1);
            d2B_dZ_dI(i,:) = d2B_dZ_dI(i,:) - sum(1. / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.,1);
            
            Ref = [1, 0 0.; 0, -1, 0; 0 0 -1]';
            for coord = 1:3
                for sc = 1:2 % sine then cos
                    for k = 1:coilData.Nt
                        idx = (i - 1) * coilData.Nt * 3 * 2 + (coord - 1) * coilData.Nt * 2 + (sc-1)*coilData.Nt + k;

                        dT_dc = zeros(coilData.M,3);
                        dxyz_dc = zeros(coilData.M,3);

                        if sc == 1
                            dxyz_dc(:,coord) = fourierData.SIN(:,k);
                            dT_dc(:,coord) = fourierData.dSIN(:,k);
                        else
                            dxyz_dc(:,coord) = fourierData.COS(:,k);
                            dT_dc(:,coord) = fourierData.dCOS(:,k);
                        end
                        dxyz_dc = dxyz_dc*Ref;
                        dT_dc = dT_dc*Ref;
                        
                        d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                        d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;
                        d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                        
                        d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, -dxyz_dc);
                        d2_dist_dR_dc = dot(repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);
                        d2_dist_dZ_dc = dot(repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);

                        term1 =  cross(dT_dc, d_diff_dR) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, d_diff_dR)./dist.^6;
                        term2 = d_distm4_dc .* dist_R .* T_x_diff + (1./dist).^4 .* (d2_dist_dR_dc .* T_x_diff + d_T_x_diff_dc .* dist_R );

                        term1z =  cross(dT_dc, d_diff_dZ) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{coilData.ts*coilData.C + i}, d_diff_dZ)./dist.^6;
                        term2z = d_distm4_dc .* dist_Z .* T_x_diff + (1./dist).^4 .* (d2_dist_dZ_dc .* T_x_diff + d_T_x_diff_dc .* dist_Z );
                        
                        dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * -coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                        d2_Bxyz_dR_dc(idx,:) = d2_Bxyz_dR_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(term1 - 3 * term2,1);
                        d2_Bxyz_dZ_dc(idx,:) = d2_Bxyz_dZ_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(term1z - 3 * term2z,1);
                    end
                end
            end
        end


        for t = 1:coilData.ts-1
            diff = point-coilData.coil_field{t*coilData.C + i};
            d_diff_dR = repmat(x_R, size(coilData.coil_tangents{t*coilData.C + i},1), 1);
            d_diff_dZ = repmat(x_Z, size(coilData.coil_tangents{t*coilData.C + i},1), 1);
                    
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
            dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
            dist_Z     = diff(:,3)./dist;

            dist3_R        = 3 * dist.^2 .* dist_R;
            dist3_Theta    = 3 * dist.^2 .* dist_Theta;
            dist3_Z        = 3 * dist.^2 .* dist_Z;
            distm4_R       =-4 * dist.^(-5) .* dist_R;
            distm4_Z       =-4 * dist.^(-5) .* dist_Z;
            d2_dist_dR2    = (dist - dist_R.* ((x/R) * diff(:,1) + (y/R) * diff(:,2))    )./(dist.^2);
            d2_dist_dR_dZ  = -diff(:,3)./dist.^2 .* dist_R;
            d2_dist_d2Z    = ( dist - dist_Z.*diff(:,3) )./ dist.^2;
            
            T_x_diff = cross(coilData.coil_tangents{t*coilData.C + i}, diff);
            T_x_diff_R = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{t*coilData.C + i},1), 1));
            T_x_diff_Theta = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{t*coilData.C + i},1), 1));
            T_x_diff_Z = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{t*coilData.C + i},1), 1));

            term1RR =  - (3./dist.^4) .* dist_R .* T_x_diff_R;
            term2RR =  distm4_R .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR2 .* T_x_diff + T_x_diff_R .* dist_R);
            d2_B_d2R = d2_B_d2R + (gamma * coilData.I(i) ) / (4. * pi) * (term1RR-3.*term2RR);

            term1RZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_R;
            term2RZ =  distm4_Z .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR_dZ .* T_x_diff + T_x_diff_Z .* dist_R);
            d2_B_dR_dZ = d2_B_dR_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (term1RZ-3.*term2RZ);
            
            term1ZZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_Z;
            term2ZZ =  distm4_Z .* dist_Z .* T_x_diff + (1./dist.^4).*(d2_dist_d2Z .* T_x_diff + T_x_diff_Z .* dist_Z);
            d2_B_d2Z = d2_B_d2Z + (gamma * coilData.I(i) ) / (4. * pi) * (term1ZZ-3.*term2ZZ);
            
            dB = dB + (gamma * coilData.I(i) ) / (4. * pi) * T_x_diff./dist.^3;
            dB_dR = dB_dR + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
            dB_dTheta = dB_dTheta + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
            dB_dZ = dB_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;
            
            dB_dI(i,:) = dB_dI(i,:) + sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
            d2B_dR_dI(i,:) = d2B_dR_dI(i,:) + sum(1. / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.,1);
            d2B_dZ_dI(i,:) = d2B_dZ_dI(i,:) + sum(1. / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.,1);
        
            Rot = [cos(t * 2 * pi / coilData.ts), - sin(t * 2 * pi / coilData.ts) 0.; sin(t * 2 * pi / coilData.ts), cos(t * 2 * pi / coilData.ts), 0; 0 0 1]';
            for coord = 1:3
                for sc = 1:2 % sine then cos
                    for k = 1:coilData.Nt
                        idx = (i - 1) * coilData.Nt * 3 * 2 + (coord - 1) * coilData.Nt * 2 + (sc-1)*coilData.Nt + k;

                        dT_dc = zeros(coilData.M,3);
                        dxyz_dc = zeros(coilData.M,3);

                        if sc == 1
                            dxyz_dc(:,coord) = fourierData.SIN(:,k);
                            dT_dc(:,coord) = fourierData.dSIN(:,k);
                        else
                            dxyz_dc(:,coord) = fourierData.COS(:,k);
                            dT_dc(:,coord) = fourierData.dCOS(:,k);
                        end
                        dxyz_dc = dxyz_dc*Rot;
                        dT_dc = dT_dc*Rot;
                        
                        d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                        d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;
                        d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                        
                        d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{t*coilData.C + i}, -dxyz_dc);
                        d2_dist_dR_dc = dot(repmat(x_R, size(coilData.coil_tangents{t*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);
                        d2_dist_dZ_dc = dot(repmat(x_Z, size(coilData.coil_tangents{t*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);

                        term1 =  cross(dT_dc, d_diff_dR) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{t*coilData.C + i}, d_diff_dR)./dist.^6;
                        term2 = d_distm4_dc .* dist_R .* T_x_diff + (1./dist).^4 .* (d2_dist_dR_dc .* T_x_diff + d_T_x_diff_dc .* dist_R );

                        term1z =  cross(dT_dc, d_diff_dZ) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{t*coilData.C + i}, d_diff_dZ)./dist.^6;
                        term2z = d_distm4_dc .* dist_Z .* T_x_diff + (1./dist).^4 .* (d2_dist_dZ_dc .* T_x_diff + d_T_x_diff_dc .* dist_Z );

                        
                        dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                        d2_Bxyz_dR_dc(idx,:) = d2_Bxyz_dR_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(term1 - 3 * term2,1);
                        d2_Bxyz_dZ_dc(idx,:) = d2_Bxyz_dZ_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(term1z - 3 * term2z,1);
                    end
                end
            end


            if coilData.ss == 1
                diff = point-coilData.coil_field{coilData.ts*coilData.C +t*coilData.C + i};
                d_diff_dR = repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1);
                d_diff_dZ = repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1);
                
                dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
                dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
                dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
                dist_Z     = diff(:,3)./dist;

                dist3_R        = 3 * dist.^2 .* dist_R;
                dist3_Theta    = 3 * dist.^2 .* dist_Theta;
                dist3_Z        = 3 * dist.^2 .* dist_Z;
                distm4_R       =-4 * dist.^(-5) .* dist_R;
                distm4_Z       =-4 * dist.^(-5) .* dist_Z;
                d2_dist_dR2    = (dist - dist_R.* ((x/R) * diff(:,1) + (y/R) * diff(:,2))    )./(dist.^2);
                d2_dist_dR_dZ  = -diff(:,3)./dist.^2 .* dist_R;
                d2_dist_d2Z    = ( dist - dist_Z.*diff(:,3) )./ dist.^2;
           
                T_x_diff = cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, diff);
                T_x_diff_R = cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1));
                T_x_diff_Theta = cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1));
                T_x_diff_Z = cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1));

                
                term1RR =  - (3./dist.^4) .* dist_R .* T_x_diff_R;
                term2RR =  distm4_R .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR2 .* T_x_diff + T_x_diff_R .* dist_R);
                d2_B_d2R = d2_B_d2R + (gamma * -coilData.I(i) ) / (4. * pi) * (term1RR-3.*term2RR);

                term1RZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_R;
                term2RZ =  distm4_Z .* dist_R .* T_x_diff + (1./dist.^4).*(d2_dist_dR_dZ .* T_x_diff + T_x_diff_Z .* dist_R);
                d2_B_dR_dZ = d2_B_dR_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (term1RZ-3.*term2RZ);
            
                term1ZZ =  - (3./dist.^4) .* dist_Z .* T_x_diff_Z;
                term2ZZ =  distm4_Z .* dist_Z .* T_x_diff + (1./dist.^4).*(d2_dist_d2Z .* T_x_diff + T_x_diff_Z .* dist_Z);
                d2_B_d2Z = d2_B_d2Z + (gamma * -coilData.I(i) ) / (4. * pi) * (term1ZZ-3.*term2ZZ);
                
                dB = dB + (gamma * -coilData.I(i) ) / (4. * pi) * T_x_diff./dist.^3;
                dB_dR = dB_dR + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
                dB_dTheta = dB_dTheta + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
                dB_dZ = dB_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;
                
                dB_dI(i,:) = dB_dI(i,:) - sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
                d2B_dR_dI(i,:) = d2B_dR_dI(i,:) - sum(1. / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.,1);
                d2B_dZ_dI(i,:) = d2B_dZ_dI(i,:) - sum(1. / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.,1);
                
                Ref = [1, 0 0.; 0, -1, 0; 0 0 -1]';
                for coord = 1:3
                    for sc = 1:2 % sine then cos
                        for k = 1:coilData.Nt
                            idx = (i - 1) * coilData.Nt * 3 * 2 + (coord - 1) * coilData.Nt * 2 + (sc-1)*coilData.Nt + k;

                            dT_dc = zeros(coilData.M,3);
                            dxyz_dc = zeros(coilData.M,3);

                            if sc == 1
                                dxyz_dc(:,coord) = fourierData.SIN(:,k);
                                dT_dc(:,coord) = fourierData.dSIN(:,k);
                            else
                                dxyz_dc(:,coord) = fourierData.COS(:,k);
                                dT_dc(:,coord) = fourierData.dCOS(:,k);
                            end
                            dxyz_dc = dxyz_dc* Rot * Ref;
                            dT_dc = dT_dc* Rot * Ref;
                            
                            d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                            d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;
                            d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                            
                            d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, -dxyz_dc);
                            d2_dist_dR_dc = dot(repmat(x_R, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);
                            d2_dist_dZ_dc = dot(repmat(x_Z, size(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i},1), 1), (-dxyz_dc .* dist - d_dist_dc .* diff)./dist.^2,2);

                            term1 =  cross(dT_dc, d_diff_dR) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, d_diff_dR)./dist.^6;
                            term2 = d_distm4_dc .* dist_R .* T_x_diff + (1./dist).^4 .* (d2_dist_dR_dc .* T_x_diff + d_T_x_diff_dc .* dist_R );

                            term1z =  cross(dT_dc, d_diff_dZ) ./ dist.^3  -   d_dist3_dc .* cross(coilData.coil_tangents{coilData.ts*coilData.C +t*coilData.C + i}, d_diff_dZ)./dist.^6;
                            term2z = d_distm4_dc .* dist_Z .* T_x_diff + (1./dist).^4 .* (d2_dist_dZ_dc .* T_x_diff + d_T_x_diff_dc .* dist_Z );

                            d2_Bxyz_dR_dc(idx,:) = d2_Bxyz_dR_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(term1 - 3 * term2,1);
                            d2_Bxyz_dZ_dc(idx,:) = d2_Bxyz_dZ_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(term1z - 3 * term2z,1);

                            
                            dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * -coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                        end
                    end
                end
            
            
            end
        end

    end
    

    
    % xyz
    dBxyz = sum(dB,1) * 2 * pi / size(dB,1);
    dBxyz_dR = sum(dB_dR,1) * 2 * pi / size(dB_dR,1);
    dBxyz_dTheta = sum(dB_dTheta,1) * 2 * pi / size(dB_dTheta,1);
    dBxyz_dZ = sum(dB_dZ,1) * 2 * pi / size(dB_dZ,1);
    
    d2_Bxyz_d2R = sum(d2_B_d2R,1) * 2 * pi / size(d2_B_d2R,1);
    d2_Bxyz_dR_dZ = sum(d2_B_dR_dZ,1) * 2 * pi / size(d2_B_dR_dZ,1);
    d2_Bxyz_d2Z = sum(d2_B_d2Z,1) * 2 * pi / size(d2_B_d2Z,1); 
    
    d2_B_RphiZ_d2R   = [cos(phi).*d2_Bxyz_d2R(1) + sin(phi).*d2_Bxyz_d2R(2) , -sin(phi).*d2_Bxyz_d2R(1) + cos(phi).*d2_Bxyz_d2R(2) , d2_Bxyz_d2R(3)];
    d2_B_RphiZ_dR_dZ = [cos(phi).*d2_Bxyz_dR_dZ(1) + sin(phi).*d2_Bxyz_dR_dZ(2) , -sin(phi).*d2_Bxyz_dR_dZ(1) + cos(phi).*d2_Bxyz_dR_dZ(2) , d2_Bxyz_dR_dZ(3)];
    d2_B_RphiZ_d2Z = [cos(phi).*d2_Bxyz_d2Z(1) + sin(phi).*d2_Bxyz_d2Z(2) , -sin(phi).*d2_Bxyz_d2Z(1) + cos(phi).*d2_Bxyz_d2Z(2) , d2_Bxyz_d2Z(3)];
    
    d2_B_RphiZ_dR_dc = [cos(phi).*d2_Bxyz_dR_dc(:,1) + sin(phi).*d2_Bxyz_dR_dc(:,2) , -sin(phi).*d2_Bxyz_dR_dc(:,1) + cos(phi).*d2_Bxyz_dR_dc(:,2) , d2_Bxyz_dR_dc(:,3)];
    d2_B_RphiZ_dZ_dc = [cos(phi).*d2_Bxyz_dZ_dc(:,1) + sin(phi).*d2_Bxyz_dZ_dc(:,2) , -sin(phi).*d2_Bxyz_dZ_dc(:,1) + cos(phi).*d2_Bxyz_dZ_dc(:,2) , d2_Bxyz_dZ_dc(:,3)];
        
    dBxyz_dI = dB_dI  * 2 * pi / size(dB,1);
    dBR_dI = (x/R) * dBxyz_dI(:,1) + (y/R) * dBxyz_dI(:,2); 
    dBphi_dI = -(y/R) * dBxyz_dI(:,1) + (x/R) * dBxyz_dI(:,2);
    dBZ_dI = dBxyz_dI(:,3);
    
    d2Bxyz_dR_dI = d2B_dR_dI * 2 * pi / size(dB,1);
    d2BR_dR_dI = (x/R) * d2Bxyz_dR_dI(:,1) + (y/R) * d2Bxyz_dR_dI(:,2); 
    d2Bphi_dR_dI = -(y/R) * d2Bxyz_dR_dI(:,1) + (x/R) * d2Bxyz_dR_dI(:,2);
    d2BZ_dR_dI = d2Bxyz_dR_dI(:,3);
    
    d2Bxyz_dZ_dI = d2B_dZ_dI * 2 * pi / size(dB,1);
    d2BR_dZ_dI = (x/R) * d2Bxyz_dZ_dI(:,1) + (y/R) * d2Bxyz_dZ_dI(:,2); 
    d2Bphi_dZ_dI = -(y/R) * d2Bxyz_dZ_dI(:,1) + (x/R) * d2Bxyz_dZ_dI(:,2);
    d2BZ_dZ_dI = d2Bxyz_dZ_dI(:,3);
    
    B = [(x/R) * dBxyz(1) + (y/R) * dBxyz(2),  -(y/R) * dBxyz(1) + (x/R) * dBxyz(2)  , dBxyz(3)];
    dB_dR = [(x/R) * dBxyz_dR(1) + (y/R) * dBxyz_dR(2),  -(y/R) * dBxyz_dR(1) + (x/R) * dBxyz_dR(2)  , dBxyz_dR(3)];
    dB_dTheta = [ (B(2) + dBxyz_dTheta(1)) * cos(phi) ,  (-B(1) + dBxyz_dTheta(2)) * sin(phi)  , dBxyz_dTheta(3)];
    dB_dZ = [(x/R) * dBxyz_dZ(1) + (y/R) * dBxyz_dZ(2),  -(y/R) * dBxyz_dZ(1) + (x/R) * dBxyz_dZ(2)  , dBxyz_dZ(3)];

    BR = B(1);
    Bphi = B(2);
    BZ = B(3);

    d_BR_d_R = dB_dR(1);
    d_Bphi_d_R = dB_dR(2);
    d_BZ_d_R = dB_dR(3);

    d_BR_d_phi = dB_dTheta(1);
    d_Bphi_d_phi = dB_dTheta(2);
    d_BZ_d_phi = dB_dTheta(3);

    d_BR_d_Z = dB_dZ(1);
    d_Bphi_d_Z = dB_dZ(2);
    d_BZ_d_Z = dB_dZ(3);
    
    dB_R_dc = dB_xyz_dc(:,1)*cos(phi)+dB_xyz_dc(:,2)*sin(phi);
    dB_phi_dc = -dB_xyz_dc(:,1)*sin(phi)+dB_xyz_dc(:,2)*cos(phi);
    dB_Z_dc = dB_xyz_dc(:,3);

    d2B_R_dc_dR = d2_B_RphiZ_dR_dc(:,1);
    d2B_phi_dc_dR = d2_B_RphiZ_dR_dc(:,2);
    d2B_Z_dc_dR = d2_B_RphiZ_dR_dc(:,3);
    
    
    d2B_R_dc_dZ = d2_B_RphiZ_dZ_dc(:,1);
    d2B_phi_dc_dZ = d2_B_RphiZ_dZ_dc(:,2);
    d2B_Z_dc_dZ = d2_B_RphiZ_dZ_dc(:,3);
    
    
    % wrt c
    d2_R_BR_Bphi_dR_dc = (1/Bphi)*dB_R_dc - (1./Bphi^2)*BR*dB_phi_dc...
                         +(R/Bphi)*d2B_R_dc_dR-(R/Bphi^2)*d_BR_d_R * dB_phi_dc...
                         -(R/Bphi^2)*d_Bphi_d_R*dB_R_dc + 2 * (R*BR/Bphi^3)*d_Bphi_d_R*dB_phi_dc ...
                         -(R*BR/Bphi^2) * d2B_phi_dc_dR;

    d2_R_BR_Bphi_dZ_dc = (R/Bphi) * d2B_R_dc_dZ - (R/Bphi^2) * d_BR_d_Z*dB_phi_dc...
                         -(R/Bphi^2) * d_Bphi_d_Z* dB_R_dc + 2 * (R/Bphi^3) * BR * d_Bphi_d_Z * dB_phi_dc...
                         -(R/Bphi^2) * BR * d2B_phi_dc_dZ;
                     
    d2_R_BZ_Bphi_dR_dc = (1/Bphi)*dB_Z_dc - (1./Bphi^2)*BZ*dB_phi_dc...
                         +(R/Bphi)*d2B_Z_dc_dR-(R/Bphi^2)*d_BZ_d_R * dB_phi_dc...
                         -(R/Bphi^2)*d_Bphi_d_R*dB_Z_dc + 2 * (R*BZ/Bphi^3)*d_Bphi_d_R*dB_phi_dc ...
                         -(R*BZ/Bphi^2) * d2B_phi_dc_dR;
                     
    d2_R_BZ_Bphi_dZ_dc = (R/Bphi) * d2B_Z_dc_dZ - (R/Bphi^2) * d_BZ_d_Z * dB_phi_dc...
                         -(R/Bphi^2) * d_Bphi_d_Z* dB_Z_dc + 2 * (R/Bphi^3) * BZ * d_Bphi_d_Z * dB_phi_dc...
                         -(R/Bphi^2) * BZ * d2B_phi_dc_dZ;
                     
    
                     
                     
                     
                     
                     
     % wrt I
    d2_R_BR_Bphi_dR_dI = (1/Bphi)*dBR_dI - (1./Bphi^2)*BR*dBphi_dI...
                         +(R/Bphi)*d2BR_dR_dI-(R/Bphi^2)*d_BR_d_R * dBphi_dI...
                         -(R/Bphi^2)*d_Bphi_d_R*dBR_dI + 2 * (R*BR/Bphi^3)*d_Bphi_d_R*dBphi_dI ...
                         -(R*BR/Bphi^2) * d2Bphi_dR_dI;

    d2_R_BR_Bphi_dZ_dI = (R/Bphi) * d2BR_dZ_dI - (R/Bphi^2) * d_BR_d_Z*dBphi_dI...
                         -(R/Bphi^2) * d_Bphi_d_Z* dBR_dI + 2 * (R/Bphi^3) * BR * d_Bphi_d_Z * dBphi_dI...
                         -(R/Bphi^2) * BR * d2Bphi_dZ_dI;
                     
    d2_R_BZ_Bphi_dR_dI = (1/Bphi)*dBZ_dI - (1./Bphi^2)*BZ*dBphi_dI...
                         +(R/Bphi)*d2BZ_dR_dI-(R/Bphi^2)*d_BZ_d_R * dBphi_dI...
                         -(R/Bphi^2)*d_Bphi_d_R*dBZ_dI + 2 * (R*BZ/Bphi^3)*d_Bphi_d_R*dBphi_dI ...
                         -(R*BZ/Bphi^2) * d2Bphi_dR_dI;
                     
    d2_R_BZ_Bphi_dZ_dI = (R/Bphi) * d2BZ_dZ_dI - (R/Bphi^2) * d_BZ_d_Z * dBphi_dI...
                         -(R/Bphi^2) * d_Bphi_d_Z* dBZ_dI + 2 * (R/Bphi^3) * BZ * d_Bphi_d_Z * dBphi_dI...
                         -(R/Bphi^2) * BZ * d2Bphi_dZ_dI;
                     
                     
      
                     
                      
                     
    d2_R_BR_Bphi_d2R =   2 * dB_dR(1) / B(2) - (2 / B(2)^2) * B(1) * dB_dR(2)...
                      + (R/B(2)) * d2_B_RphiZ_d2R(1) - 2 * (R / B(2)^2)*dB_dR(1)*dB_dR(2) ...
                      + 2*(R/B(2)^3)*B(1) * dB_dR(2)^2 - (R/B(2)^2)*B(1) * d2_B_RphiZ_d2R(2);

    d2_R_BR_Bphi_dZ_dR = dB_dZ(1) / B(2) - (1/B(2)^2) * B(1) * dB_dZ(2) ...
                         + (R/B(2)) * d2_B_RphiZ_dR_dZ(1) - (R/B(2)^2) * dB_dR(1) * dB_dZ(2) - (R/B(2)^2) * dB_dZ(1) * dB_dR(2)...
                         + 2 * (R/B(2)^3) * B(1) * dB_dR(2) * dB_dZ(2) - (R / B(2)^2) * B(1) * d2_B_RphiZ_dR_dZ(2);
                     
    d2_R_BZ_Bphi_d2R =    2 * dB_dR(3) / B(2) - (2 / B(2)^2) * B(3) * dB_dR(2)...
                      + (R/B(2)) * d2_B_RphiZ_d2R(3) - 2 * (R / B(2)^2)*dB_dR(3)*dB_dR(2) ...
                      + 2*(R/B(2)^3)*B(3) * dB_dR(2)^2 - (R/B(2)^2)*B(3) * d2_B_RphiZ_d2R(2);
                     
    d2_R_BZ_Bphi_dZ_dR = dB_dZ(3) / B(2) - (1/B(2)^2) * B(3) * dB_dZ(2) ...
                         + (R/B(2)) * d2_B_RphiZ_dR_dZ(3) - (R/B(2)^2) * dB_dR(3) * dB_dZ(2) - (R/B(2)^2) * dB_dZ(3) * dB_dR(2)...
                         + 2 * (R/B(2)^3) * B(3) * dB_dR(2) * dB_dZ(2) - (R / B(2)^2) * B(3) * d2_B_RphiZ_dR_dZ(2);
                     
                     
    d2_R_BR_Bphi_d2Z = (R / B(2) ) * d2_B_RphiZ_d2Z(1) - 2 * ( R / B(2)^2 ) *  dB_dZ(1) * dB_dZ(2)...
                      + 2*(R/B(2)^3)*B(1) * dB_dZ(2)^2 - (R / B(2)^2) * B(1) * d2_B_RphiZ_d2Z(2);
    d2_R_BZ_Bphi_d2Z = (R / B(2) ) * d2_B_RphiZ_d2Z(3) - 2 * ( R / B(2)^2 ) *  dB_dZ(3) * dB_dZ(2)...
                      + 2*(R/B(2)^3)*B(3) * dB_dZ(2)^2 - (R / B(2)^2) * B(3) * d2_B_RphiZ_d2Z(2);
end

