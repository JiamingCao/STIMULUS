function [factors] = GetFactorFromResistorNet(num_x, num_y, num_z, Ymat, Az, anode_loc, cathode_loc, direction) % polar coordinates
% Solver of problem Ymat*V=I using lsqr
% electrode locations are matrices, each row is one electrode [theta, phi]
dx = 2/(num_x - 1);
dy = 2/(num_y - 1);
dz = 2/(num_z - 1);
D = [dx, dy, dz]';
num_nodes = num_x*num_y*num_z;
factors = zeros(num_nodes, size(anode_loc, 1));
[x1, y1, z1] = polar_to_cartesian(1, anode_loc(:,1), anode_loc(:,2));
[x2, y2, z2] = polar_to_cartesian(1, cathode_loc(:,1), cathode_loc(:,2));
x1_idx = ceil(num_x/2) + sign(x1).*floor((abs(x1)+eps)/dx);
y1_idx = ceil(num_y/2) + sign(y1).*floor((abs(y1)+eps)/dy);
z1_idx = ceil(num_z/2) + sign(z1).*floor((abs(z1)+eps)/dz);

x2_idx = ceil(num_x/2) + sign(x2).*floor((abs(x2)+eps)/dx);
y2_idx = ceil(num_y/2) + sign(y2).*floor((abs(y2)+eps)/dy);
z2_idx = ceil(num_z/2) + sign(z2).*floor((abs(z2)+eps)/dz);

anode_ind = sub2ind([num_x, num_y, num_z], x1_idx, y1_idx, z1_idx);
cathode_ind = sub2ind([num_x, num_y, num_z], x2_idx, y2_idx, z2_idx);
Imat = zeros(num_nodes, length(anode_ind));
for i = 1:length(anode_ind)
    I = zeros(num_nodes, 1);
    I(cathode_ind(i)) = -1; % cathode position
    I(anode_ind(i)) = 1; % anode position
    Imat(:,i) = I;
end
d = D(direction);
try
    Ymat = gpuArray(Ymat);
    Imat = gpuArray(Imat);
    Az = gpuArray(Az);
catch ;
end
parfor i = 1:length(anode_ind)
    V_at_node = lsqr(Ymat, Imat(:,i), [], 50000);
    Vzz = (Az * V_at_node);
    factors(:, i) = gather((Vzz)/d^2);
end
