clear;
tic;
%% Calculate frequencies
n_pairs = 2;% fc=2.5k; 2pairs df=15, 415; 4 pairs df=10, 290; 8 pairs df=10, 160; 16 pairs df=10, 77-optimal

f_center = 2000;      % center frequency
delta_f = 40;
if n_pairs == 2
%     freq_shift = [0;40];
    freq_shift = [0;50];
else
%     freq_shift = [0, 10, 23, 10, 0, -10, -23, -10]';
%     freq_shift = [0, 5, 10, 15, 30, 15, 10, 5, 0, -5, -10, -15, -30, -15, -10, -5]';
%     freq_shift = 0.25*delta_f * (1 - cos(4*pi/n_pairs * ([1:n_pairs/2] - 1)));
%     freq_shift = [freq_shift, -freq_shift]';
    freq_shift = [0:2*delta_f/n_pairs:delta_f/2,delta_f/2-2*delta_f/n_pairs:-2*delta_f/n_pairs:-delta_f/2,-delta_f/2+2*delta_f/n_pairs:2*delta_f/n_pairs:-2*delta_f/n_pairs]';
end
freq = f_center + freq_shift;
%% Discretize space and calculate factors
CenterAmplitude = 440; % change this (uA/cm^2)
Fs = 25;
vstart = -70;
twidth = 1000;
% sources
phi_step = 2*pi/n_pairs;
phi = (0:phi_step:2*pi-phi_step)';
theta1 = 1*pi/6;
theta2 = pi - theta1;
ElectrodeLoc_upper = [theta1*ones(size(phi)), phi];
ElectrodeLoc_lower = [theta2*ones(size(phi)), phi];

num_x = 51;
num_y = 51;
num_z = 51;
xgrid = linspace(-1,1,num_x); 
ygrid = linspace(-1,1,num_y);
zgrid = linspace(-1,1,num_z);
[X, Y, Z] = meshgrid(xgrid, ygrid, zgrid);
radius1 = 0.86; % brain/skull/scalp ~= 7.9/8.6/9.2
radius2 = 0.93;
radius3 = 1;
unit_len_equals = 1; %cm
R = [1, 15, 1]'; % Change skull resistance here
num_nodes = length(xgrid)*length(ygrid)*length(zgrid);

nodes = [X(:), Y(:), Z(:)];
nodes_in_sphere = 3*(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius3^2+eps);
nodes_in_sphere(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius2^2+eps) = 2;
nodes_in_sphere(nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <= radius1^2+eps) = 1;

lst_add = find(nodes_in_sphere > 0);
lst_add = [lst_add, lst_add+1;lst_add, lst_add+length(xgrid); lst_add, lst_add+length(xgrid)*length(ygrid)];
lst_add(lst_add(:,2)>num_nodes, :) = [];
lst_add(nodes_in_sphere(lst_add(:,2))<1, :) = [];
seperation = find(diff(lst_add(:,1))<0);
lst_add_x = [lst_add(1:seperation(1), :), ones(seperation(1), 1)];
lst_add_y = [lst_add(seperation(1)+1:seperation(2),:), ones(seperation(2)-seperation(1), 1)];
lst_add_z = [lst_add(seperation(2)+1:end,:), ones(length(lst_add)-seperation(2), 1)];
lst_add = [lst_add, -1./R(max(nodes_in_sphere(lst_add), [], 2))];

Ax = sparse(lst_add_x(:,1), lst_add_x(:,2), lst_add_x(:,3), num_nodes, num_nodes);
Ax = Ax + Ax';
Ax = Ax - 2*speye(size(Ax));

Ay = sparse(lst_add_y(:,1), lst_add_y(:,2), lst_add_y(:,3), num_nodes, num_nodes);
Ay = Ay + Ay';
Ay = Ay - 2*speye(size(Ay));

Az = sparse(lst_add_z(:,1), lst_add_z(:,2), lst_add_z(:,3), num_nodes, num_nodes);
Az = Az + Az';
Az = Az - 2*speye(size(Az));

Ymat = sparse(lst_add(:,1), lst_add(:,2), lst_add(:,3), num_nodes, num_nodes);
Ymat = Ymat + Ymat'; % make the matrix symmetric
Ymat = Ymat - diag(sum(Ymat)); % and deal with the diagonals

clear lst* % save some memory

% Get factors, and scale them properly. "unit_len_equals" is not useful
% here.
factors = GetFactorFromResistorNet(num_x, num_y, num_z, Ymat, Az, ElectrodeLoc_upper, ElectrodeLoc_lower, 3);
focus = 0.4;
[~, idx_focus] = min(abs(zgrid - focus));
factors = factors./factors(sub2ind(size(X), ceil(num_x/2),ceil(num_y/2),idx_focus),:)*CenterAmplitude;
% factors = factors * scale;
% currentmag = currentmag/currentmag(sub2ind(size(X),ceil(num_x/2),ceil(num_y/2),ceil(num_z/2)),1)*CenterAmplitude;
roi = nodes_in_sphere == 1;% & nodes(:,3)>0.36 & nodes(:,3)<0.44;
factors2 = factors(roi,:);
fire = nan(size(factors,1), 1);
fire2 = nan(size(factors2, 1), 1);
[b,a]=butter(3,[100/(Fs*500), 1000/(Fs*500)]);
thresh = 30;
parfor idx = 1:size(factors2,1)
    amp = factors2(idx,:)';
%     [T,S] = HodgkinHuxleyPV_fast(vstart,twidth,freq,amp,Fs);
%     [T,S] = HodgkinHuxleyE_iso(vstart,twidth,freq,amp,Fs);
    [T,S] = HodgkinHuxleyFast(vstart,twidth,freq,amp,Fs);
%     [T,S] = Sejnowski(vstart,twidth,freq,amp,Fs);
    if length(S) < twidth*Fs+1
           fire2(idx) = -1;
        continue;
    end

    Vfilt = filtfilt(b,a,S);
    if max(Vfilt(10000:20000))>thresh
        fire2(idx) = 1;
    else
        fire2(idx) = 0;
    end        
end
fire(roi) = fire2;
fire = reshape(fire, num_x, num_y, num_z);
% fire = permute(fire, [2,1,3]);
toc;

% figure, imagesc(xgrid, ygrid, squeeze(fire(:,:,36))');
figure, slice(xgrid,ygrid,zgrid,permute(fire, [2,1,3]),0,0,0.4);

% save('Results/TI_2pair_15_01')
