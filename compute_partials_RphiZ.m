function [BR,Bphi,BZ,...
          d_BR_d_R,d_Bphi_d_R,d_BZ_d_R,...
          d_BR_d_phi,d_Bphi_d_phi,d_BZ_d_phi,...
          d_BR_d_Z,d_Bphi_d_Z,d_BZ_d_Z] = compute_partials_RphiZ(R,phi,Z, coilData)
    
    x = R * cos(phi);
    y = R * sin(phi);
    z = Z;

    
    x_R = [x/R, y/R, 0];
    x_Theta= [-y, x, 0];
    x_Z = [0, 0, 1];
    

    %% Declaration of variables
    gamma = 4 * pi * 10^(-7) ;
    dB = zeros(coilData.M,3);
    dB_dR = zeros(coilData.M,3);
    dB_dTheta = zeros(coilData.M,3);
    dB_dZ = zeros(coilData.M,3);
    point = [x y z];


    for i = 1: coilData.C
        %% Numerical integration of Biot-Savart law
        diff = point-coilData.coil_field{i};
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
        dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
        dist_Z     = diff(:,3)./dist;

        dist3_R        = 3 * dist.^2 .* dist_R;
        dist3_Theta    = 3 * dist.^2 .* dist_Theta;
        dist3_Z        = 3 * dist.^2 .* dist_Z;

        T_x_diff = cross(coilData.coil_tangents{i}, diff);
        T_x_diff_R = cross(coilData.coil_tangents{i}, repmat(x_R, size(coilData.coil_tangents{i},1), 1));
        T_x_diff_Theta = cross(coilData.coil_tangents{i}, repmat(x_Theta, size(coilData.coil_tangents{i},1), 1));
        T_x_diff_Z = cross(coilData.coil_tangents{i}, repmat(x_Z, size(coilData.coil_tangents{i},1), 1));
        
        dB = dB + (gamma * coilData.I(i) ) / (4. * pi) * cross(coilData.coil_tangents{i}, diff)./dist.^3;
        dB_dR = dB_dR + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
        dB_dTheta = dB_dTheta + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
        dB_dZ = dB_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;


        if coilData.ss == 1
            diff = point-coilData.coil_field{coilData.nfp*coilData.C + i};
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
            dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
            dist_Z     = diff(:,3)./dist;

            dist3_R        = 3 * dist.^2 .* dist_R;
            dist3_Theta    = 3 * dist.^2 .* dist_Theta;
            dist3_Z        = 3 * dist.^2 .* dist_Z;

            T_x_diff = cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, diff);
            T_x_diff_R = cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{coilData.nfp*coilData.C + i},1), 1));
            T_x_diff_Theta = cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{coilData.nfp*coilData.C + i},1), 1));
            T_x_diff_Z = cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{coilData.nfp*coilData.C + i},1), 1));
            
            dB = dB + (gamma * -coilData.I(i) ) / (4. * pi) * cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, diff)./dist.^3;
            dB_dR = dB_dR + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
            dB_dTheta = dB_dTheta + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
            dB_dZ = dB_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;
        end








        for t = 1:coilData.nfp-1
            diff = point-coilData.coil_field{t*coilData.C + i};
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
            dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
            dist_Z     = diff(:,3)./dist;

            dist3_R        = 3 * dist.^2 .* dist_R;
            dist3_Theta    = 3 * dist.^2 .* dist_Theta;
            dist3_Z        = 3 * dist.^2 .* dist_Z;

            T_x_diff = cross(coilData.coil_tangents{t*coilData.C + i}, diff);
            T_x_diff_R = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{t*coilData.C + i},1), 1));
            T_x_diff_Theta = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{t*coilData.C + i},1), 1));
            T_x_diff_Z = cross(coilData.coil_tangents{t*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{t*coilData.C + i},1), 1));
            
            dB = dB + (gamma * coilData.I(i)) / (4. * pi) * cross(coilData.coil_tangents{t*coilData.C + i}, diff)./dist.^3;
            dB_dR = dB_dR + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
            dB_dTheta = dB_dTheta + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
            dB_dZ = dB_dZ + (gamma * coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;



            if coilData.ss == 1
                diff = point-coilData.coil_field{coilData.nfp*coilData.C +t*coilData.C + i};
                dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
                dist_R     = ((x/R) * diff(:,1) + (y/R) * diff(:,2) )./dist;
                dist_Theta = (-y * diff(:,1) + x * diff(:,2))./dist;
                dist_Z     = diff(:,3)./dist;

                dist3_R        = 3 * dist.^2 .* dist_R;
                dist3_Theta    = 3 * dist.^2 .* dist_Theta;
                dist3_Z        = 3 * dist.^2 .* dist_Z;

                T_x_diff = cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, diff);
                T_x_diff_R = cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, repmat(x_R, size(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},1), 1));
                T_x_diff_Theta = cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, repmat(x_Theta, size(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},1), 1));
                T_x_diff_Z = cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, repmat(x_Z, size(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},1), 1));
                
                dB = dB + (gamma * -coilData.I(i)) / (4. * pi) * cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, diff)./dist.^3;
                dB_dR = dB_dR + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_R .* dist.^3 - dist3_R.*T_x_diff)./dist.^6.;
                dB_dTheta = dB_dTheta + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Theta .* dist.^3 - dist3_Theta.*T_x_diff)./dist.^6.;
                dB_dZ = dB_dZ + (gamma * -coilData.I(i) ) / (4. * pi) * (T_x_diff_Z .* dist.^3 - dist3_Z.*T_x_diff)./dist.^6.;
            end
        end

    end
    % xyz
    dBxyz = sum(dB,1) * 2 * pi / size(dB,1);
    dBxyz_dR = sum(dB_dR,1) * 2 * pi / size(dB_dR,1);
    dBxyz_dTheta = sum(dB_dTheta,1) * 2 * pi / size(dB_dTheta,1);
    dBxyz_dZ = sum(dB_dZ,1) * 2 * pi / size(dB_dZ,1);
    
    
    
    
    
    
    
    
    
    
    B = [(x/R) * dBxyz(1) + (y/R) * dBxyz(2),  -(y/R) * dBxyz(1) + (x/R) * dBxyz(2)  , dBxyz(3)];
    BR = B(1);
    Bphi = B(2);
    BZ = B(3);
    dB_dR = [(x/R) * dBxyz_dR(1) + (y/R) * dBxyz_dR(2),  -(y/R) * dBxyz_dR(1) + (x/R) * dBxyz_dR(2)  , dBxyz_dR(3)];
    dB_dTheta = [ (B(2) + dBxyz_dTheta(1)) * cos(phi) ,  (-B(1) + dBxyz_dTheta(2)) * sin(phi)  , dBxyz_dTheta(3)];
    dB_dZ = [(x/R) * dBxyz_dZ(1) + (y/R) * dBxyz_dZ(2),  -(y/R) * dBxyz_dZ(1) + (x/R) * dBxyz_dZ(2)  , dBxyz_dZ(3)];

    d_BR_d_R = dB_dR(1);
    d_Bphi_d_R = dB_dR(2);
    d_BZ_d_R = dB_dR(3);

    d_BR_d_phi = dB_dTheta(1);
    d_Bphi_d_phi = dB_dTheta(2);
    d_BZ_d_phi = dB_dTheta(3);

    d_BR_d_Z = dB_dZ(1);
    d_Bphi_d_Z = dB_dZ(2);
    d_BZ_d_Z = dB_dZ(3);
end