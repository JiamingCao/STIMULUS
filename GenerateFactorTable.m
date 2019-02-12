clear;
tic;
num_x = 51; % Define discretization level here
num_y = 51;
num_z = 51;
x_mid = ceil(num_x/2);
z_mid = ceil(num_z/2);
y_mid = ceil(num_y/2);

xgrid = linspace(-1,1,num_x); 
ygrid = linspace(-1,1,num_y);
zgrid = linspace(-1,1,num_z);

[X, Y, Z] = meshgrid(xgrid, ygrid, zgrid);
radius1 = 0.86; % brain/skull/scalp ~= 7.9/8.6/9.2
radius2 = 0.93;
radius3 = 1;
unit_len_equals = 1; %cm
R = [1, 15, 1]'; % Resistance of three layers, i.e. sigma_1, sigma_2, and sigma_3
num_nodes = length(xgrid)*length(ygrid)*length(zgrid);

% lst_add = 3*zeros(num_nodes, 3);
% lst_add_z = 3*zeros(num_nodes, 3);
% lst_add_x = 3*zeros(num_nodes, 3);
% lst_add_y = 3*zeros(num_nodes, 3);
% cnt = 1;
% cnt_x = 1;
% cnt_y = 1;
% cnt_z = 1;

nodes = [X(:), Y(:), Z(:)];
nodes_in_sphere = 3*(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius3^2+eps);
nodes_in_sphere(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius2^2+eps) = 2;
nodes_in_sphere(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius1^2+eps) = 1;

% for i = 1:num_nodes
%     if nodes_in_sphere(i)
%         if i+1 <= num_nodes && nodes_in_sphere(i+1) % (i+1,j,k)
%             lst_add(cnt,:) = [i, i+1, -1/R(max(nodes_in_sphere(i),nodes_in_sphere(i+1)))];
%             lst_add_x(cnt_x,:) = [i, i+1, 1];
%             cnt=cnt+1;
%             cnt_x = cnt_x+1;
%         end
%         if i+length(xgrid) <= num_nodes && nodes_in_sphere(i+length(xgrid)) % (i,j+1,k)
%              lst_add(cnt,:) = [i, i+length(xgrid), -1/R(max(nodes_in_sphere(i),nodes_in_sphere(i+length(xgrid))))];
%              lst_add_y(cnt_y,:) = [i, i+length(xgrid), 1];
%              cnt=cnt+1;
%              cnt_y=cnt_y+1;
%         end
%         if i+length(xgrid)*length(ygrid) <= num_nodes && nodes_in_sphere(i+length(xgrid)*length(ygrid)) % (i,j,k+1)
%             lst_add(cnt,:) = [i, i+length(xgrid)*length(ygrid), -1/R(max(nodes_in_sphere(i),nodes_in_sphere(i+length(xgrid)*length(ygrid))))];
%             lst_add_z(cnt_z,:) = [i, i+length(xgrid)*length(ygrid), 1];
%             cnt=cnt+1;
%             cnt_z = cnt_z + 1;
%         end
%     end
% end   

lst_add = find(nodes_in_sphere > 0);
lst_add = [lst_add, lst_add+1;lst_add, lst_add+length(xgrid); lst_add, lst_add+length(xgrid)*length(ygrid)];
lst_add(lst_add(:,2)>num_nodes, :) = [];    % get rid of false connections to nodes outside of the cube
lst_add(nodes_in_sphere(lst_add(:,2))<1, :) = [];   % get rid of false connections to nodes inside of the cube but outside of the brain
seperation = find(diff(lst_add(:,1))<0);    % get the separation idx for x-, y-, z- directional connections
lst_add_x = [lst_add(1:seperation(1), :), ones(seperation(1), 1)];
lst_add_y = [lst_add(seperation(1)+1:seperation(2),:), ones(seperation(2)-seperation(1), 1)];
lst_add_z = [lst_add(seperation(2)+1:end,:), ones(length(lst_add)-seperation(2), 1)];
lst_add = [lst_add, -1./R(max(nodes_in_sphere(lst_add), [], 2))];   % set up the Kirchhoff equations

Ax = sparse(lst_add_x(:,1), lst_add_x(:,2), lst_add_x(:,3), num_nodes, num_nodes);
Ax = Ax + Ax';
Ax = Ax - 2*speye(size(Ax));

Ay = sparse(lst_add_y(:,1), lst_add_y(:,2), lst_add_y(:,3), num_nodes, num_nodes);
Ay = Ay + Ay';
Ay = Ay - 2*speye(size(Ay));

Az = sparse(lst_add_z(:,1), lst_add_z(:,2), lst_add_z(:,3), num_nodes, num_nodes);
Az = Az + Az';
Az = Az - 2*speye(size(Az));    % second derivative: the factor of next and previous node is 1, factor of itself is -2

Ymat = sparse(lst_add(:,1), lst_add(:,2), lst_add(:,3), num_nodes, num_nodes);
Ymat = Ymat + Ymat'; % make the matrix symmetric
Ymat = Ymat - diag(sum(Ymat)); % and deal with the diagonals

clear lst* % save some memory

toc;
locations = cell(x_mid * z_mid, 1); % [i,j,k]
cnt = 1;

dy1 = ygrid(2) - ygrid(1);
for z_idx = z_mid+1:num_z
    r_square = radius3^2 - zgrid(z_idx)^2 + eps;
    for x_idx = x_mid:num_x
        y_idx = y_mid + floor((sqrt(r_square - xgrid(x_idx)^2)+eps)/dy1);
        if x_idx > y_idx
            break;
        end
        locations{cnt} = num2str([x_idx, y_idx, z_idx]);
        cnt = cnt + 1;
        if x_idx == y_idx
            break;
        end
    end
end
locations = locations(~cellfun('isempty',locations));
factors_xx = cell(size(locations));
factors_yy = cell(size(locations));
factors_zz = cell(size(locations));

dx = unit_len_equals*2/(num_x - 1);
dy = unit_len_equals*2/(num_y - 1);
dz = unit_len_equals*2/(num_z - 1);
for i = 1:length(locations)
    tmp = str2num(locations{i});
    x = tmp(1);
    y = tmp(2);
    z = tmp(3);
    anode_ind = sub2ind(size(X), x, y, z);
    cathode_ind = sub2ind(size(X), x, y, num_z+1-z);
    I = zeros(num_nodes, 1);
    I(anode_ind) = 1;
    I(cathode_ind) = -1;
    Imat(:,i) = I;
end
try
    Ymat = gpuArray(Ymat);
    Imat = gpuArray(Imat);
    Az = gpuArray(Az);
catch;
end
parfor i = 1:length(locations)
    V_at_node = lsqr(Ymat, Imat(:,i), [], 50000);
    Vzz = (Az * V_at_node);
    Vxx = Ax * V_at_node;
    Vyy = Ay * V_at_node;
    factors_xx{i} = gather((Vxx)/dx^2);
    factors_yy{i} = gather((Vyy)/dy^2);
    factors_zz{i} = gather((Vzz)/dz^2);
end

FactorTable = containers.Map(locations, factors_zz);
save('FactorTable_51_15zz2', 'FactorTable');
toc
