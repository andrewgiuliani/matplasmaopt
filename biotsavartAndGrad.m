%% Biot-Savart integration on a generic curve
function [B1, B2, B3,...
          dB1_dx, dB2_dx, dB3_dx,...
          dB1_dy, dB2_dy, dB3_dy,...
          dB1_dz, dB2_dz, dB3_dz]=biotsavartAndGrad(x,y,z, coilData)

coil_points = coilData.coil_points;
I = coilData.I;
coil_tangents = coilData.coil_tangents;
nCoils = coilData.C;
nfp = coilData.nfp;
ss = coilData.ss;


%% Declaration of variables
gamma = 4 * pi * 10^(-7) ;


dB = zeros(size(coil_points,1),3);
dB_dx = zeros(size(coil_points,1),3);
dB_dy = zeros(size(coil_points,1),3);
dB_dz = zeros(size(coil_points,1),3);
point = [x y z];



for i = 1: nCoils
    %% Numerical integration of Biot-Savart law
    
    diff = point-coilData.coil_field{i};
    dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
    
    dist_x = diff(:,1)./dist;
    dist_y = diff(:,2)./dist;
    dist_z = diff(:,3)./dist;
    
    d_point_dx = repmat([1,0,0], size(coilData.coil_tangents{i},1), 1);
    d_point_dy = repmat([0,1,0], size(coilData.coil_tangents{i},1), 1);
    d_point_dz = repmat([0,0,1], size(coilData.coil_tangents{i},1), 1);
    
    dB = dB + (gamma * I(i) ) / (4. * pi) * cross(coil_tangents{i}, diff)./dist.^3;
    dB_dx = dB_dx + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dx).*dist.^3 ...
                                                   -3*dist.^2 .* dist_x .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
    dB_dy = dB_dy + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dy).*dist.^3 ...
                                                   -3*dist.^2 .* dist_y .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
    dB_dz = dB_dz + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{i}, d_point_dz).*dist.^3 ...
                                                   -3*dist.^2 .* dist_z .* cross(coil_tangents{i}, diff) ) ./ dist.^6;
    if ss == 1
        diff = point-coilData.coil_field{coilData.nfp*coilData.C + i};
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        
        dist_x = diff(:,1)./dist;
        dist_y = diff(:,2)./dist;
        dist_z = diff(:,3)./dist;
        
        dB = dB + (gamma * -I(i) ) / (4. * pi) * cross(coil_tangents{nfp*nCoils + i}, diff)./dist.^3;
        dB_dx = dB_dx + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dx).*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
        dB_dy = dB_dy + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dy).*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
        dB_dz = dB_dz + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils + i}, d_point_dz).*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* cross(coil_tangents{nfp*nCoils + i}, diff) ) ./ dist.^6;
    end 
    
    for t = 1:nfp-1
        diff = point-coilData.coil_field{t*coilData.C + i};
        dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
        
        dist_x = diff(:,1)./dist;
        dist_y = diff(:,2)./dist;
        dist_z = diff(:,3)./dist;
        
        dB = dB + (gamma * I(i)) / (4. * pi) * cross(coil_tangents{t*nCoils + i}, diff)./dist.^3;
        dB_dx = dB_dx + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dx).*dist.^3 ...
                                                       -3*dist.^2 .* dist_x .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
        dB_dy = dB_dy + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dy).*dist.^3 ...
                                                       -3*dist.^2 .* dist_y .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
        dB_dz = dB_dz + (gamma * I(i) ) / (4. * pi) * ( cross(coil_tangents{t*nCoils + i}, d_point_dz).*dist.^3 ...
                                                       -3*dist.^2 .* dist_z .* cross(coil_tangents{t*nCoils + i}, diff) ) ./ dist.^6;
                                                
        if ss == 1
            diff = point-coilData.coil_field{coilData.nfp*coilData.C +t*coilData.C + i};
            dist = sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            
            dist_x = diff(:,1)./dist;
            dist_y = diff(:,2)./dist;
            dist_z = diff(:,3)./dist;
            
            dB = dB + (gamma * -I(i)) / (4. * pi) * cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff)./dist.^3;
            dB_dx = dB_dx + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dx).*dist.^3 ...
                                                           -3*dist.^2 .* dist_x .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
            dB_dy = dB_dy + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dy).*dist.^3 ...
                                                           -3*dist.^2 .* dist_y .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
            dB_dz = dB_dz + (gamma * -I(i) ) / (4. * pi) * ( cross(coil_tangents{nfp*nCoils +t*nCoils + i}, d_point_dz).*dist.^3 ...
                                                           -3*dist.^2 .* dist_z .* cross(coil_tangents{nfp*nCoils +t*nCoils + i}, diff) ) ./ dist.^6;
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

end
