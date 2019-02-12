function [ElectrodeAmplitudeVector, currentmag, cancel_mag, exitflag] = spatial_interfere(focus, cancel, factors, MaxVol, num_x, num_y, num_z) % Cartesian coordinates
x_mid = ceil(num_x/2);
y_mid = ceil(num_y/2);
z_mid = ceil(num_z/2);

dx = 2/(num_x - 1);
dy = 2/(num_y - 1);
dz = 2/(num_z - 1);

focus_ind = sub2ind([num_x, num_y, num_z], round(focus(:,1)/dx)+x_mid, round(focus(:,2)/dy)+y_mid, round(focus(:,3)/dz)+z_mid);
cancel_ind = sub2ind([num_x, num_y, num_z], round(cancel(:,1)/dx)+x_mid, round(cancel(:,2)/dy)+y_mid, round(cancel(:,3)/dz)+z_mid);
A_mat = factors(cancel_ind, :)';
a0 = factors(focus_ind, :)';
Q = A_mat * A_mat';

% Q = Q + 0.0025* (factors'*factors)/(num_x*num_y*num_z);

options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','MaxIterations',30000, 'Display', 'off');
Desired_focus_mag =ones(1,size(focus,1));
num_pairs = size(factors, 2);
tempvec = ones(num_pairs,1);% MaxVol = 0.1;
[ElectrodeAmplitudeVector,~,exitflag] = quadprog(Q,[],[],[],a0', Desired_focus_mag, -MaxVol*tempvec,MaxVol*tempvec,[],options);%-30*tempvec,30*tempvec
exitflag = ~(exitflag==1);
if exitflag
    ElectrodeAmplitudeVector = zeros(num_pairs, 1);
end
cancel_mag = abs(ElectrodeAmplitudeVector'*A_mat);
% attained_focus_mag = a0'*ElectrodeAmplitudeVector;
currentmag = abs(factors*ElectrodeAmplitudeVector);

