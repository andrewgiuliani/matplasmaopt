%% Biot-Savart integration on a generic curve
function [B1, B2, B3,...
          dB1_dx, dB1_dy, dB1_dz,...
          dB2_dx, dB2_dy, dB2_dz,...
          dB3_dx, dB3_dy, dB3_dz,...
          dB1_dc, dB2_dc, dB3_dc,...
          dB1_dI, dB2_dI, dB3_dI,...
          dB1x_dI,dB1y_dI,dB1z_dI,...
          dB2x_dI,dB2y_dI,dB2z_dI,...
          dB3x_dI,dB3y_dI,dB3z_dI,...
          dB1_dRcos, dB2_dRcos, dB3_dRcos,...
          dB1_dZsin, dB2_dZsin, dB3_dZsin,...
          dB1x_dRcos,dB1y_dRcos,dB1z_dRcos,...
          dB2x_dRcos,dB2y_dRcos,dB2z_dRcos,...
          dB3x_dRcos,dB3y_dRcos,dB3z_dRcos,...
          dB1x_dZsin,dB1y_dZsin,dB1z_dZsin,...
          dB2x_dZsin,dB2y_dZsin,dB2z_dZsin,...
          dB3x_dZsin,dB3y_dZsin,dB3z_dZsin,...
          B1x_dc, B1y_dc, B1z_dc, ...
          B2x_dc, B2y_dc, B2z_dc, ...
          B3x_dc, B3y_dc, B3z_dc]=biotsavartwithpartials(x,y,z, coilData, ...
                                                                  x_Rcos, y_Rcos, z_Rcos,...
                                                                  x_Zsin, y_Zsin, z_Zsin)

coil_points = coilData.coil_points;
I = coilData.I;
coil_tangents = coilData.coil_tangents;
nCoils = coilData.C;
nfp = coilData.nfp;
ss = coilData.ss;
fourierData = coilData.fourierData;

%% Declaration of variables
gamma = 4 * pi * 10^(-7) ;


dB = zeros(size(coil_points,1),3);
dB_dx = zeros(size(coil_points,1),3);
dB_dy = zeros(size(coil_points,1),3);
dB_dz = zeros(size(coil_points,1),3);
dB_xyz_dc = zeros(numel(coilData.coil_coeffs),3);
dB_dI = zeros(size(coilData.I,1),3);
dBx_dI = zeros(size(coilData.I,1),3);
dBy_dI = zeros(size(coilData.I,1),3);
dBz_dI = zeros(size(coilData.I,1),3);
dB_xyz_dRcos = zeros(numel(x_Rcos),3);
dB_xyz_dZsin = zeros(numel(x_Zsin),3);

d2_Bxyz_coils_dx_dc = zeros(numel(coilData.coil_coeffs),3);
d2_Bxyz_coils_dy_dc = zeros(numel(coilData.coil_coeffs),3);
d2_Bxyz_coils_dz_dc = zeros(numel(coilData.coil_coeffs),3);

dB_dx_dRcos = zeros(numel(x_Rcos),3);
dB_dx_dZsin = zeros(numel(x_Zsin),3);

dB_dy_dRcos = zeros(numel(x_Rcos),3);
dB_dy_dZsin = zeros(numel(x_Zsin),3);

dB_dz_dRcos = zeros(numel(x_Rcos),3);
dB_dz_dZsin = zeros(numel(x_Zsin),3);

point = [x y z];
    
d_point_dx = repmat([1,0,0], size(coilData.coil_tangents{1},1), 1);
d_point_dy = repmat([0,1,0], size(coilData.coil_tangents{1},1), 1);
d_point_dz = repmat([0,0,1], size(coilData.coil_tangents{1},1), 1);

for i = 1: nCoils
    %% Numerical integration of Biot-Savart law
    
    diff = point-coilData.coil_field{i};
    dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
    
    dist_x = diff(:,1)./dist;
    dist_y = diff(:,2)./dist;
    dist_z = diff(:,3)./dist;

    
    T_x_diff = cross(coilData.coil_tangents{i}, diff);
    T_x_d_point_dx = cross(coil_tangents{i}, d_point_dx);
    T_x_d_point_dy = cross(coil_tangents{i}, d_point_dy);
    T_x_d_point_dz = cross(coil_tangents{i}, d_point_dz);
    
    dB = dB + (gamma * I(i) ) / (4. * pi) * cross(coil_tangents{i}, diff)./dist.^3;
    dB_dx = dB_dx + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dx).*dist.^3 ...
                                                   -3*dist.^2 .* dist_x .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
    dB_dy = dB_dy + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dy).*dist.^3 ...
                                                   -3*dist.^2 .* dist_y .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
    dB_dz = dB_dz + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dz).*dist.^3 ...
                                                   -3*dist.^2 .* dist_z .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
                                               
    dB_dI(i,:) = dB_dI(i,:) + sum(1. / (4. * pi) * T_x_diff./dist.^3, 1);                                           
    dBx_dI(i,:) = dBx_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dx.*dist.^3 ...
                                                   -3*dist.^2 .* dist_x .* T_x_diff ) ./ dist.^6, 1);  
    dBy_dI(i,:) = dBy_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dy.*dist.^3 ...
                                                   -3*dist.^2 .* dist_y .* T_x_diff ) ./ dist.^6, 1);  
    dBz_dI(i,:) = dBz_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dz.*dist.^3 ...
                                                   -3*dist.^2 .* dist_z .* T_x_diff ) ./ dist.^6, 1);
    
    

    
    
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
                d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{i}, -dxyz_dc);
                
                
                d2_dist_dxyz_dc = (dist .* (-dxyz_dc) - d_dist_dc .* diff) ./ dist.^2;
                
                
                termA_dx = (cross(dT_dc, d_point_dx) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{i},d_point_dx))./dist.^6;
                termA_dy = (cross(dT_dc, d_point_dy) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{i},d_point_dy))./dist.^6;
                termA_dz = (cross(dT_dc, d_point_dz) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{i},d_point_dz))./dist.^6;
                
                termB_dx = d_distm4_dc .* T_x_diff .* dist_x + (1./dist).^4 .* (d2_dist_dxyz_dc(:,1) .* T_x_diff + d_T_x_diff_dc .* dist_x );
                termB_dy = d_distm4_dc .* T_x_diff .* dist_y + (1./dist).^4 .* (d2_dist_dxyz_dc(:,2) .* T_x_diff + d_T_x_diff_dc .* dist_y );
                termB_dz = d_distm4_dc .* T_x_diff .* dist_z + (1./dist).^4 .* (d2_dist_dxyz_dc(:,3) .* T_x_diff + d_T_x_diff_dc .* dist_z );
                
                
                d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{i}, -dxyz_dc);
                dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
            
                d2_Bxyz_coils_dx_dc(idx,:) = d2_Bxyz_coils_dx_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dx - 3 * termB_dx,1);
                d2_Bxyz_coils_dy_dc(idx,:) = d2_Bxyz_coils_dy_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dy - 3 * termB_dy,1);
                d2_Bxyz_coils_dz_dc(idx,:) = d2_Bxyz_coils_dz_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dz - 3 * termB_dz,1);
                
                
            end
        end
    end   
    
    
    

    

    
    

    for idx = 1:numel(x_Rcos)
        d_point_dRcos = repmat([x_Rcos(idx), y_Rcos(idx), z_Rcos(idx)], coilData.M, 1);
        dist_Rcos     = (diff(:,1)*x_Rcos(idx) + diff(:,2)*y_Rcos(idx) + diff(:,3)*z_Rcos(idx) )./dist;

        
        d_point_dZsin = repmat([x_Zsin(idx), y_Zsin(idx), z_Zsin(idx)], coilData.M, 1);
        dist_Zsin     = (diff(:,1)*x_Zsin(idx) + diff(:,2)*y_Zsin(idx) + diff(:,3)*z_Zsin(idx) )./dist;


        dB_xyz_dRcos(idx,:) = dB_xyz_dRcos(idx,:) + (2 * pi / coilData.M) * sum( (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dRcos).*dist.^3 ...
                                                                                         -3*dist.^2 .* dist_Rcos .* cross(coil_tangents{i}, diff) ) ./ dist.^6, 1);
        dB_xyz_dZsin(idx,:) = dB_xyz_dZsin(idx,:) + (2 * pi / coilData.M) * sum( (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dZsin).*dist.^3 ...
                                                                                         -3*dist.^2 .* dist_Zsin .* cross(coil_tangents{i}, diff) ) ./ dist.^6, 1);
    
    

                                                                                     
                                                                                     
                                                                                     
                                                                                     
                                                                                     
                                                                                     
        T_x_d_point_dRcos = cross(coil_tangents{i}, d_point_dRcos);     
        T_x_d_point_dZsin = cross(coil_tangents{i}, d_point_dZsin); 
                                                                                     
        dist_Rcos = dot(d_point_dRcos, diff, 2)./dist;
        dist_Zsin = dot(d_point_dZsin, diff, 2)./dist;

        dist_x_Rcos = (x_Rcos(idx) .* dist - dist_Rcos.*diff(:,1))./dist.^2;
        dist_x_Zsin = (x_Zsin(idx) .* dist - dist_Zsin.*diff(:,1))./dist.^2;

        dist_y_Rcos = (y_Rcos(idx) .* dist - dist_Rcos.*diff(:,2))./dist.^2;
        dist_y_Zsin = (y_Zsin(idx) .* dist - dist_Zsin.*diff(:,2))./dist.^2;

        dist_z_Rcos = (z_Rcos(idx) .* dist - dist_Rcos.*diff(:,3))./dist.^2;
        dist_z_Zsin = (z_Zsin(idx) .* dist - dist_Zsin.*diff(:,3))./dist.^2;

        term1x_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Rcos;
        term1x_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Zsin;

        term1y_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Rcos;
        term1y_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Zsin;

        term1z_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Rcos;
        term1z_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Zsin; 
        
        term2x_dRcos  =   (dist_x_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dRcos;
        term2x_dZsin  =   (dist_x_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dZsin; 
        
        term2y_dRcos  =   (dist_y_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dRcos;
        term2y_dZsin  =   (dist_y_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dZsin; 
        
        term2z_dRcos  =   (dist_z_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dRcos;
        term2z_dZsin  =   (dist_z_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dZsin; 
        
        
        dB_dx_dRcos(idx,:) = dB_dx_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1x_dRcos -3 * term2x_dRcos, 1);
        dB_dx_dZsin(idx,:) = dB_dx_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1x_dZsin -3 * term2x_dZsin, 1);
    
        dB_dy_dRcos(idx,:) = dB_dy_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1y_dRcos -3 * term2y_dRcos, 1);
        dB_dy_dZsin(idx,:) = dB_dy_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1y_dZsin -3 * term2y_dZsin, 1);

        dB_dz_dRcos(idx,:) = dB_dz_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1z_dRcos -3 * term2z_dRcos, 1);
        dB_dz_dZsin(idx,:) = dB_dz_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1z_dZsin -3 * term2z_dZsin, 1);
    end
                                               
    if ss == 1
        diff = point-coilData.coil_field{coilData.nfp*coilData.C + i};
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        
        dist_x = diff(:,1)./dist;
        dist_y = diff(:,2)./dist;
        dist_z = diff(:,3)./dist;
        
        T_x_diff = cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, diff);
        T_x_d_point_dx = cross(coil_tangents{nfp*nCoils + i}, d_point_dx);
        T_x_d_point_dy = cross(coil_tangents{nfp*nCoils + i}, d_point_dy);
        T_x_d_point_dz = cross(coil_tangents{nfp*nCoils + i}, d_point_dz);
        
        dB = dB + (gamma * -I(i) ) / (4. * pi) * cross(coil_tangents{nfp*nCoils + i}, diff)./dist.^3;
        dB_dx = dB_dx + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dx).*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
        dB_dy = dB_dy + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dy).*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
        dB_dz = dB_dz + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dz).*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
                                                   
        dB_dI(i,:) = dB_dI(i,:) - sum(1. / (4. * pi) * T_x_diff./dist.^3,1);                                          
        dBx_dI(i,:) = dBx_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dx.*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* T_x_diff ) ./ dist.^6, 1);  
        dBy_dI(i,:) = dBy_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dy.*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* T_x_diff ) ./ dist.^6, 1);  
        dBz_dI(i,:) = dBz_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dz.*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* T_x_diff ) ./ dist.^6, 1);     
                                                   
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

                    d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;
                    d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{coilData.nfp*coilData.C + i}, -dxyz_dc);

                    dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * -coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                
                
                    
                    
                    d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                    d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                    d2_dist_dxyz_dc = (dist .* (-dxyz_dc) - d_dist_dc .* diff) ./ dist.^2;
                
                
                    termA_dx = (cross(dT_dc, d_point_dx) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C + i},d_point_dx))./dist.^6;
                    termA_dy = (cross(dT_dc, d_point_dy) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C + i},d_point_dy))./dist.^6;
                    termA_dz = (cross(dT_dc, d_point_dz) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C + i},d_point_dz))./dist.^6;

                    termB_dx = d_distm4_dc .* T_x_diff .* dist_x + (1./dist).^4 .* (d2_dist_dxyz_dc(:,1) .* T_x_diff + d_T_x_diff_dc .* dist_x );
                    termB_dy = d_distm4_dc .* T_x_diff .* dist_y + (1./dist).^4 .* (d2_dist_dxyz_dc(:,2) .* T_x_diff + d_T_x_diff_dc .* dist_y );
                    termB_dz = d_distm4_dc .* T_x_diff .* dist_z + (1./dist).^4 .* (d2_dist_dxyz_dc(:,3) .* T_x_diff + d_T_x_diff_dc .* dist_z );

                    d2_Bxyz_coils_dx_dc(idx,:) = d2_Bxyz_coils_dx_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dx - 3 * termB_dx,1);
                    d2_Bxyz_coils_dy_dc(idx,:) = d2_Bxyz_coils_dy_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dy - 3 * termB_dy,1);
                    d2_Bxyz_coils_dz_dc(idx,:) = d2_Bxyz_coils_dz_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dz - 3 * termB_dz,1);
                
                
                
                
                
                end
            end
        end                               

        

              

        for idx = 1:numel(x_Rcos)
            d_point_dRcos = repmat([x_Rcos(idx), y_Rcos(idx), z_Rcos(idx)], coilData.M, 1);
            dist_Rcos     = (diff(:,1)*x_Rcos(idx) + diff(:,2)*y_Rcos(idx) + diff(:,3)*z_Rcos(idx) )./dist;

             

            d_point_dZsin = repmat([x_Zsin(idx), y_Zsin(idx), z_Zsin(idx)], coilData.M, 1);
            dist_Zsin     = (diff(:,1)*x_Zsin(idx) + diff(:,2)*y_Zsin(idx) + diff(:,3)*z_Zsin(idx) )./dist;


            dB_xyz_dRcos(idx,:) = dB_xyz_dRcos(idx,:) + (2 * pi / coilData.M) * sum( (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dRcos).*dist.^3 ...
                                                                                     -3*dist.^2 .* dist_Rcos .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6, 1);
            dB_xyz_dZsin(idx,:) = dB_xyz_dZsin(idx,:) + (2 * pi / coilData.M) * sum( (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dZsin).*dist.^3 ...
                                                                                     -3*dist.^2 .* dist_Zsin .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6, 1);
        
        
        
        
            T_x_d_point_dRcos = cross(coil_tangents{nfp*nCoils + i}, d_point_dRcos);     
            T_x_d_point_dZsin = cross(coil_tangents{nfp*nCoils + i}, d_point_dZsin); 

            dist_Rcos = dot(d_point_dRcos, diff, 2)./dist;
            dist_Zsin = dot(d_point_dZsin, diff, 2)./dist;

            dist_x_Rcos = (x_Rcos(idx) .* dist - dist_Rcos.*diff(:,1))./dist.^2;
            dist_x_Zsin = (x_Zsin(idx) .* dist - dist_Zsin.*diff(:,1))./dist.^2;

            dist_y_Rcos = (y_Rcos(idx) .* dist - dist_Rcos.*diff(:,2))./dist.^2;
            dist_y_Zsin = (y_Zsin(idx) .* dist - dist_Zsin.*diff(:,2))./dist.^2;

            dist_z_Rcos = (z_Rcos(idx) .* dist - dist_Rcos.*diff(:,3))./dist.^2;
            dist_z_Zsin = (z_Zsin(idx) .* dist - dist_Zsin.*diff(:,3))./dist.^2;

            term1x_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Rcos;
            term1x_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Zsin;

            term1y_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Rcos;
            term1y_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Zsin;

            term1z_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Rcos;
            term1z_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Zsin; 

            term2x_dRcos  =   (dist_x_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dRcos;
            term2x_dZsin  =   (dist_x_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dZsin; 

            term2y_dRcos  =   (dist_y_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dRcos;
            term2y_dZsin  =   (dist_y_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dZsin; 

            term2z_dRcos  =   (dist_z_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dRcos;
            term2z_dZsin  =   (dist_z_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dZsin; 


            dB_dx_dRcos(idx,:) = dB_dx_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1x_dRcos -3 * term2x_dRcos, 1);
            dB_dx_dZsin(idx,:) = dB_dx_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1x_dZsin -3 * term2x_dZsin, 1);

            dB_dy_dRcos(idx,:) = dB_dy_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1y_dRcos -3 * term2y_dRcos, 1);
            dB_dy_dZsin(idx,:) = dB_dy_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1y_dZsin -3 * term2y_dZsin, 1);

            dB_dz_dRcos(idx,:) = dB_dz_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1z_dRcos -3 * term2z_dRcos, 1);
            dB_dz_dZsin(idx,:) = dB_dz_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1z_dZsin -3 * term2z_dZsin, 1);
        
        
        
        end

        
        
    end 
    
    for t = 1:nfp-1
        diff = point-coilData.coil_field{t*coilData.C + i};
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        
        dist_x = diff(:,1)./dist;
        dist_y = diff(:,2)./dist;
        dist_z = diff(:,3)./dist;
        
        T_x_diff = cross(coilData.coil_tangents{t*coilData.C + i}, diff);
        T_x_d_point_dx = cross(coil_tangents{t*nCoils + i}, d_point_dx);
        T_x_d_point_dy = cross(coil_tangents{t*nCoils + i}, d_point_dy);
        T_x_d_point_dz = cross(coil_tangents{t*nCoils + i}, d_point_dz);
        
        dB = dB + (gamma * I(i)) / (4. * pi) * cross(coil_tangents{t*nCoils + i}, diff)./dist.^3;
        dB_dx = dB_dx + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dx).*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
        dB_dy = dB_dy + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dy).*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
        dB_dz = dB_dz + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dz).*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
        dB_dI(i,:) = dB_dI(i,:) + sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
        dBx_dI(i,:) = dBx_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dx.*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* T_x_diff ) ./ dist.^6, 1);  
        dBy_dI(i,:) = dBy_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dy.*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* T_x_diff ) ./ dist.^6, 1);  
        dBz_dI(i,:) = dBz_dI(i,:) + sum(1. / (4. * pi) * ( T_x_d_point_dz.*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* T_x_diff ) ./ dist.^6, 1);  
                                                   
        Rot = [cos(t * 2 * pi / coilData.nfp), - sin(t * 2 * pi / coilData.nfp) 0.; sin(t * 2 * pi / coilData.nfp), cos(t * 2 * pi / coilData.nfp), 0; 0 0 1]';
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
                    d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;

                    d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{t*coilData.C + i}, -dxyz_dc);
                    dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                    
                    
                    
                    d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                    d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                    d2_dist_dxyz_dc = (dist .* (-dxyz_dc) - d_dist_dc .* diff) ./ dist.^2;
                
                
                    termA_dx = (cross(dT_dc, d_point_dx) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{t*coilData.C + i},d_point_dx))./dist.^6;
                    termA_dy = (cross(dT_dc, d_point_dy) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{t*coilData.C + i},d_point_dy))./dist.^6;
                    termA_dz = (cross(dT_dc, d_point_dz) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{t*coilData.C + i},d_point_dz))./dist.^6;

                    termB_dx = d_distm4_dc .* T_x_diff .* dist_x + (1./dist).^4 .* (d2_dist_dxyz_dc(:,1) .* T_x_diff + d_T_x_diff_dc .* dist_x );
                    termB_dy = d_distm4_dc .* T_x_diff .* dist_y + (1./dist).^4 .* (d2_dist_dxyz_dc(:,2) .* T_x_diff + d_T_x_diff_dc .* dist_y );
                    termB_dz = d_distm4_dc .* T_x_diff .* dist_z + (1./dist).^4 .* (d2_dist_dxyz_dc(:,3) .* T_x_diff + d_T_x_diff_dc .* dist_z );

                    d2_Bxyz_coils_dx_dc(idx,:) = d2_Bxyz_coils_dx_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dx - 3 * termB_dx,1);
                    d2_Bxyz_coils_dy_dc(idx,:) = d2_Bxyz_coils_dy_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dy - 3 * termB_dy,1);
                    d2_Bxyz_coils_dz_dc(idx,:) = d2_Bxyz_coils_dz_dc(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum(termA_dz - 3 * termB_dz,1);
                    
                    
                    
                end
            end
        end
        
        

        
        
        
        
        for idx = 1:numel(x_Rcos)
            d_point_dRcos = repmat([x_Rcos(idx), y_Rcos(idx), z_Rcos(idx)], coilData.M, 1);
            dist_Rcos     = (diff(:,1)*x_Rcos(idx) + diff(:,2)*y_Rcos(idx) + diff(:,3)*z_Rcos(idx) )./dist;


            
            d_point_dZsin = repmat([x_Zsin(idx), y_Zsin(idx), z_Zsin(idx)], coilData.M, 1);
            dist_Zsin     = (diff(:,1)*x_Zsin(idx) + diff(:,2)*y_Zsin(idx) + diff(:,3)*z_Zsin(idx) )./dist;


            dB_xyz_dRcos(idx,:) = dB_xyz_dRcos(idx,:) + (2 * pi / coilData.M) * sum( (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dRcos).*dist.^3 ...
                                                                                     -3*dist.^2 .* dist_Rcos .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6, 1);
            dB_xyz_dZsin(idx,:) = dB_xyz_dZsin(idx,:) + (2 * pi / coilData.M) * sum( (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dZsin).*dist.^3 ...
                                                                                     -3*dist.^2 .* dist_Zsin .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6, 1);
        
        
        
        
        
        
        
            T_x_d_point_dRcos = cross(coil_tangents{t*nCoils + i}, d_point_dRcos);     
            T_x_d_point_dZsin = cross(coil_tangents{t*nCoils + i}, d_point_dZsin); 

            dist_Rcos = dot(d_point_dRcos, diff, 2)./dist;
            dist_Zsin = dot(d_point_dZsin, diff, 2)./dist;

            dist_x_Rcos = (x_Rcos(idx) .* dist - dist_Rcos.*diff(:,1))./dist.^2;
            dist_x_Zsin = (x_Zsin(idx) .* dist - dist_Zsin.*diff(:,1))./dist.^2;

            dist_y_Rcos = (y_Rcos(idx) .* dist - dist_Rcos.*diff(:,2))./dist.^2;
            dist_y_Zsin = (y_Zsin(idx) .* dist - dist_Zsin.*diff(:,2))./dist.^2;

            dist_z_Rcos = (z_Rcos(idx) .* dist - dist_Rcos.*diff(:,3))./dist.^2;
            dist_z_Zsin = (z_Zsin(idx) .* dist - dist_Zsin.*diff(:,3))./dist.^2;

            term1x_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Rcos;
            term1x_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Zsin;

            term1y_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Rcos;
            term1y_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Zsin;

            term1z_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Rcos;
            term1z_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Zsin; 

            term2x_dRcos  =   (dist_x_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dRcos;
            term2x_dZsin  =   (dist_x_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dZsin; 

            term2y_dRcos  =   (dist_y_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dRcos;
            term2y_dZsin  =   (dist_y_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dZsin; 

            term2z_dRcos  =   (dist_z_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dRcos;
            term2z_dZsin  =   (dist_z_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dZsin; 


            dB_dx_dRcos(idx,:) = dB_dx_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1x_dRcos -3 * term2x_dRcos, 1);
            dB_dx_dZsin(idx,:) = dB_dx_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1x_dZsin -3 * term2x_dZsin, 1);

            dB_dy_dRcos(idx,:) = dB_dy_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1y_dRcos -3 * term2y_dRcos, 1);
            dB_dy_dZsin(idx,:) = dB_dy_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1y_dZsin -3 * term2y_dZsin, 1);

            dB_dz_dRcos(idx,:) = dB_dz_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1z_dRcos -3 * term2z_dRcos, 1);
            dB_dz_dZsin(idx,:) = dB_dz_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (coilData.I(i) / (4 * pi) ) * sum( term1z_dZsin -3 * term2z_dZsin, 1);
       end
        
                                                   
        if ss == 1
            diff = point-coilData.coil_field{coilData.nfp*coilData.C +t*coilData.C + i};
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            
            dist_x = diff(:,1)./dist;
            dist_y = diff(:,2)./dist;
            dist_z = diff(:,3)./dist;
            
            T_x_diff = cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, diff);
            T_x_d_point_dx = cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dx);
            T_x_d_point_dy = cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dy);
            T_x_d_point_dz = cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dz); 
            
            dB = dB + (gamma * -I(i)) / (4. * pi) * cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff)./dist.^3;
            dB_dx = dB_dx + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dx).*dist.^3 ...
                                                           -3*dist.^2 .* dist_x .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
            dB_dy = dB_dy + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dy).*dist.^3 ...
                                                           -3*dist.^2 .* dist_y .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
            dB_dz = dB_dz + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dz).*dist.^3 ...
                                                           -3*dist.^2 .* dist_z .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
        
            dB_dI(i,:) = dB_dI(i,:) - sum(1. / (4. * pi) * T_x_diff./dist.^3,1);
            dBx_dI(i,:) = dBx_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dx.*dist.^3 ...
                                                           -3*dist.^2 .* dist_x .* T_x_diff ) ./ dist.^6, 1);  
            dBy_dI(i,:) = dBy_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dy.*dist.^3 ...
                                                           -3*dist.^2 .* dist_y .* T_x_diff ) ./ dist.^6, 1);  
            dBz_dI(i,:) = dBz_dI(i,:) - sum(1. / (4. * pi) * ( T_x_d_point_dz.*dist.^3 ...
                                                           -3*dist.^2 .* dist_z .* T_x_diff ) ./ dist.^6, 1);
                                                       
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

                        d_dist3_dc = 3*dist .* dot(diff, -dxyz_dc, 2) ;

                        d_T_x_diff_dc = cross(dT_dc, diff) + cross(coilData.coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i}, -dxyz_dc);
                        dB_xyz_dc(idx,:) = dB_xyz_dc(idx,:) + (2 * pi / coilData.M) * sum((gamma * -coilData.I(i) ) / (4. * pi) *  ( d_T_x_diff_dc.*dist.^3 - d_dist3_dc.*T_x_diff ) ./dist.^6, 1);
                        
                        
                        
                        d_dist_dc  = dot(diff, -dxyz_dc, 2) ./ dist;
                        d_distm4_dc = -4*dist.^(-6) .* dot(diff, -dxyz_dc, 2) ;
                        d2_dist_dxyz_dc = (dist .* (-dxyz_dc) - d_dist_dc .* diff) ./ dist.^2;


                        termA_dx = (cross(dT_dc, d_point_dx) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},d_point_dx))./dist.^6;
                        termA_dy = (cross(dT_dc, d_point_dy) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},d_point_dy))./dist.^6;
                        termA_dz = (cross(dT_dc, d_point_dz) .* dist.^3 - d_dist3_dc .* cross(coil_tangents{coilData.nfp*coilData.C +t*coilData.C + i},d_point_dz))./dist.^6;

                        termB_dx = d_distm4_dc .* T_x_diff .* dist_x + (1./dist).^4 .* (d2_dist_dxyz_dc(:,1) .* T_x_diff + d_T_x_diff_dc .* dist_x );
                        termB_dy = d_distm4_dc .* T_x_diff .* dist_y + (1./dist).^4 .* (d2_dist_dxyz_dc(:,2) .* T_x_diff + d_T_x_diff_dc .* dist_y );
                        termB_dz = d_distm4_dc .* T_x_diff .* dist_z + (1./dist).^4 .* (d2_dist_dxyz_dc(:,3) .* T_x_diff + d_T_x_diff_dc .* dist_z );

                        d2_Bxyz_coils_dx_dc(idx,:) = d2_Bxyz_coils_dx_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dx - 3 * termB_dx,1);
                        d2_Bxyz_coils_dy_dc(idx,:) = d2_Bxyz_coils_dy_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dy - 3 * termB_dy,1);
                        d2_Bxyz_coils_dz_dc(idx,:) = d2_Bxyz_coils_dz_dc(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum(termA_dz - 3 * termB_dz,1);
                        
                        
                    end
                end
            end
        

            
          
            
            
            
            for idx = 1:numel(x_Rcos)
                d_point_dRcos = repmat([x_Rcos(idx), y_Rcos(idx), z_Rcos(idx)], coilData.M, 1);
                dist_Rcos     = (diff(:,1)*x_Rcos(idx) + diff(:,2)*y_Rcos(idx) + diff(:,3)*z_Rcos(idx) )./dist;

                
                
                d_point_dZsin = repmat([x_Zsin(idx), y_Zsin(idx), z_Zsin(idx)], coilData.M, 1);
                dist_Zsin     = (diff(:,1)*x_Zsin(idx) + diff(:,2)*y_Zsin(idx) + diff(:,3)*z_Zsin(idx) )./dist;


                dB_xyz_dRcos(idx,:) = dB_xyz_dRcos(idx,:) + (2 * pi / coilData.M) * sum( (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dRcos).*dist.^3 ...
                                                                                          -3*dist.^2 .* dist_Rcos .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6, 1);
                dB_xyz_dZsin(idx,:) = dB_xyz_dZsin(idx,:) + (2 * pi / coilData.M) * sum( (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dZsin).*dist.^3 ...
                                                                                          -3*dist.^2 .* dist_Zsin .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6, 1);
            
            
            
            
            
            
            
                T_x_d_point_dRcos = cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dRcos);     
                T_x_d_point_dZsin = cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dZsin); 

                dist_Rcos = dot(d_point_dRcos, diff, 2)./dist;
                dist_Zsin = dot(d_point_dZsin, diff, 2)./dist;

                dist_x_Rcos = (x_Rcos(idx) .* dist - dist_Rcos.*diff(:,1))./dist.^2;
                dist_x_Zsin = (x_Zsin(idx) .* dist - dist_Zsin.*diff(:,1))./dist.^2;

                dist_y_Rcos = (y_Rcos(idx) .* dist - dist_Rcos.*diff(:,2))./dist.^2;
                dist_y_Zsin = (y_Zsin(idx) .* dist - dist_Zsin.*diff(:,2))./dist.^2;

                dist_z_Rcos = (z_Rcos(idx) .* dist - dist_Rcos.*diff(:,3))./dist.^2;
                dist_z_Zsin = (z_Zsin(idx) .* dist - dist_Zsin.*diff(:,3))./dist.^2;

                term1x_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Rcos;
                term1x_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dx .* dist_Zsin;

                term1y_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Rcos;
                term1y_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dy .* dist_Zsin;

                term1z_dRcos =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Rcos;
                term1z_dZsin =   -3.*(1./dist).^4 .* T_x_d_point_dz .* dist_Zsin; 

                term2x_dRcos  =   (dist_x_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dRcos;
                term2x_dZsin  =   (dist_x_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_x./dist.^5) .* T_x_diff + (dist_x./dist.^4) .* T_x_d_point_dZsin; 

                term2y_dRcos  =   (dist_y_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dRcos;
                term2y_dZsin  =   (dist_y_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_y./dist.^5) .* T_x_diff + (dist_y./dist.^4) .* T_x_d_point_dZsin; 

                term2z_dRcos  =   (dist_z_Rcos ./ dist.^4 - 4 * dist_Rcos .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dRcos;
                term2z_dZsin  =   (dist_z_Zsin ./ dist.^4 - 4 * dist_Zsin .* dist_z./dist.^5) .* T_x_diff + (dist_z./dist.^4) .* T_x_d_point_dZsin; 


                dB_dx_dRcos(idx,:) = dB_dx_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1x_dRcos -3 * term2x_dRcos, 1);
                dB_dx_dZsin(idx,:) = dB_dx_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1x_dZsin -3 * term2x_dZsin, 1);

                dB_dy_dRcos(idx,:) = dB_dy_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1y_dRcos -3 * term2y_dRcos, 1);
                dB_dy_dZsin(idx,:) = dB_dy_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1y_dZsin -3 * term2y_dZsin, 1);

                dB_dz_dRcos(idx,:) = dB_dz_dRcos(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1z_dRcos -3 * term2z_dRcos, 1);
                dB_dz_dZsin(idx,:) = dB_dz_dZsin(idx,:) + (2 * pi / coilData.M) * gamma * (-coilData.I(i) / (4 * pi) ) * sum( term1z_dZsin -3 * term2z_dZsin, 1);
        
        
        
            
            
            end
            
            
        
        end
    end 
        
end
Bfield = sum(dB,1) * 2 * pi / size(dB,1);
B1 = Bfield(1);
B2 = Bfield(2);
B3 = Bfield(3);

B_dx = sum(dB_dx,1) * 2 * pi / size(dB_dx,1);
B_dy = sum(dB_dy,1) * 2 * pi / size(dB_dy,1);
B_dz = sum(dB_dz,1) * 2 * pi / size(dB_dz,1);

dB1_dx = B_dx(1);
dB1_dy = B_dy(1);
dB1_dz = B_dz(1);

dB2_dx = B_dx(2);
dB2_dy = B_dy(2);
dB2_dz = B_dz(2);

dB3_dx = B_dx(3);
dB3_dy = B_dy(3);
dB3_dz = B_dz(3);

dB1_dc = dB_xyz_dc(:,1);
dB2_dc = dB_xyz_dc(:,2);
dB3_dc = dB_xyz_dc(:,3);

dB1_dI = dB_dI(:,1) * 2 * pi / size(dB,1);
dB2_dI = dB_dI(:,2) * 2 * pi / size(dB,1);
dB3_dI = dB_dI(:,3) * 2 * pi / size(dB,1);

dB1x_dI = dBx_dI(:,1)* 2 * pi / size(dB,1);
dB1y_dI = dBy_dI(:,1)* 2 * pi / size(dB,1);
dB1z_dI = dBz_dI(:,1)* 2 * pi / size(dB,1);

dB2x_dI = dBx_dI(:,2)* 2 * pi / size(dB,1);
dB2y_dI = dBy_dI(:,2)* 2 * pi / size(dB,1);
dB2z_dI = dBz_dI(:,2)* 2 * pi / size(dB,1);

dB3x_dI = dBx_dI(:,3)* 2 * pi / size(dB,1);
dB3y_dI = dBy_dI(:,3)* 2 * pi / size(dB,1);
dB3z_dI = dBz_dI(:,3)* 2 * pi / size(dB,1);

dB1_dRcos = dB_xyz_dRcos(:,1);
dB1_dZsin = dB_xyz_dZsin(:,1);

dB2_dRcos = dB_xyz_dRcos(:,2);
dB2_dZsin = dB_xyz_dZsin(:,2);

dB3_dRcos = dB_xyz_dRcos(:,3);
dB3_dZsin = dB_xyz_dZsin(:,3);

B1x_dc = d2_Bxyz_coils_dx_dc(:,1);
B2x_dc = d2_Bxyz_coils_dx_dc(:,2);
B3x_dc = d2_Bxyz_coils_dx_dc(:,3);

B1y_dc = d2_Bxyz_coils_dy_dc(:,1);
B2y_dc = d2_Bxyz_coils_dy_dc(:,2);
B3y_dc = d2_Bxyz_coils_dy_dc(:,3);

B1z_dc = d2_Bxyz_coils_dz_dc(:,1);
B2z_dc = d2_Bxyz_coils_dz_dc(:,2);
B3z_dc = d2_Bxyz_coils_dz_dc(:,3);

dB1x_dRcos = dB_dx_dRcos(:,1);
dB1y_dRcos = dB_dy_dRcos(:,1);
dB1z_dRcos = dB_dz_dRcos(:,1);
dB2x_dRcos = dB_dx_dRcos(:,2);
dB2y_dRcos = dB_dy_dRcos(:,2);
dB2z_dRcos = dB_dz_dRcos(:,2);
dB3x_dRcos = dB_dx_dRcos(:,3);
dB3y_dRcos = dB_dy_dRcos(:,3);
dB3z_dRcos = dB_dz_dRcos(:,3);


dB1x_dZsin = dB_dx_dZsin(:,1);
dB1y_dZsin = dB_dy_dZsin(:,1);
dB1z_dZsin = dB_dz_dZsin(:,1);
dB2x_dZsin = dB_dx_dZsin(:,2);
dB2y_dZsin = dB_dy_dZsin(:,2);
dB2z_dZsin = dB_dz_dZsin(:,2);
dB3x_dZsin = dB_dx_dZsin(:,3);
dB3y_dZsin = dB_dy_dZsin(:,3);
dB3z_dZsin = dB_dz_dZsin(:,3);
end
