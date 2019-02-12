tic;
clear;
rng('shuffle');
% load('FactorTable_51_15zz.mat');
% if isempty(gcp('nocreate'))
%     myCluster = parcluster('local');
%     core = myCluster.NumWorkers;
%     parpool(core);
% end
for take=1:5
    load('FactorTable_51_80zz2.mat');
    num_x = 51;
    num_y = 51;
    num_z = 51;
    x_mid = ceil(num_x/2);
    y_mid = ceil(num_y/2);
    z_mid = ceil(num_z/2);
    num_nodes = num_x*num_y*num_z;
    xgrid = linspace(-1,1,num_x); 
    ygrid = linspace(-1,1,num_y);
    zgrid = linspace(-1,1,num_z);
    dx = xgrid(2) - xgrid(1);
    dy = ygrid(2) - ygrid(1);
    dz = zgrid(2) - zgrid(1);
    [X, Y, Z] = meshgrid(xgrid, ygrid, zgrid);
    r_brain = 0.86; % brain/skull/scalp ~= 7.9/8.6/9.2
    nodes = [X(:), Y(:), Z(:)];

    locations = keys(FactorTable);
    locations = str2num(cell2mat(locations')); % 1/8 of a circle
    tmp = locations;
    tmp(:,1:2) = [num_x, num_y] - locations(:,1:2) +1;
    surface = unique([locations;tmp], 'rows'); % get the diagonal

    surface = unique([surface; [surface(:,2), surface(:,1), surface(:,3)]], 'rows'); % diagonal quadrants
    tmp = [surface(:,1), 1+num_y - surface(:,2), surface(:,3)];
    surface = unique([surface; tmp], 'rows'); % the other quadrants, and finish

    x_surf = xgrid(surface(:,1));
    y_surf = ygrid(surface(:,2));
    z_surf = zgrid(surface(:,3));
    r = sqrt(x_surf.^2+y_surf.^2+z_surf.^2);
    theta = acos(z_surf ./ r);
    phi = atan2(y_surf, x_surf);

    % Change cancel points and focus points here
    n_patches = 8;
    theta_patch{1} = [pi/12, pi/2];
    phi_patch{1} = [-pi/8, pi/8];
    num_pairs{1} = 300;

    theta_patch{2} = [pi/12, pi/2];
    phi_patch{2} = [pi/4-pi/8 pi/4+pi/8];
    num_pairs{2} = 300;
    
    theta_patch{3} = [pi/12, pi/2];
    phi_patch{3} = [pi/2-pi/8, pi/2+pi/8];
    num_pairs{3} = 300;
    
    theta_patch{4} = [pi/12, pi/2];
    phi_patch{4} = [3*pi/4-pi/8, 3*pi/4+pi/8];
    num_pairs{4} = 300;
    
    theta_patch{5} = [pi/12, pi/2];
    phi_patch{5} = [pi-pi/8, pi+pi/8];
    num_pairs{5} = 300;
    
    theta_patch{6} = [pi/12, pi/2];
    phi_patch{6} = [5*pi/4-pi/8, 5*pi/4+pi/8];
    num_pairs{6} = 300;
    
    theta_patch{7} = [pi/12, pi/2];
    phi_patch{7} = [3*pi/2-pi/8, 3*pi/2+pi/8];
    num_pairs{7} = 300;
    
    theta_patch{8} = [pi/12, pi/2];
    phi_patch{8} = [7*pi/4-pi/8, 7*pi/4+pi/8];
    num_pairs{8} = 300;

    focus1{1} = [0,0,0.4]; % Cartesian
    cancel1{1} = [0,0,0.6;]; % Cartesian
    max_vol1{1} = 1.6;

    focus1{2} = [0,0,0.4]; % Cartesian
    cancel1{2} = [0,0,0.6;]; % Cartesian
    max_vol1{2} = 1.6;

    focus1{3} = [0,0,0.4]; % Cartesian
    cancel1{3} = [0,0,0.6;]; % Cartesian
    max_vol1{3} = 1.6;

    focus1{4} = [0,0,0.4]; % Cartesian
    cancel1{4} = [0,0,0.6;]; % Cartesian
    max_vol1{4} = 1.6;

    focus1{5} = [0,0,0.4]; % Cartesian
    cancel1{5} = [0,0,0.6;]; % Cartesian
    max_vol1{5} = 1.6;

    focus1{6} = [0,0,0.4]; % Cartesian
    cancel1{6} = [0,0,0.6;]; % Cartesian
    max_vol1{6} = 1.6;

    focus1{7} = [0,0,0.4]; % Cartesian
    cancel1{7} = [0,0,0.6;]; % Cartesian
    max_vol1{7} = 1.6;

    focus1{8} = [0,0,0.4]; % Cartesian
    cancel1{8} = [0,0,0.6;]; % Cartesian
    max_vol1{8} = 1.6;

    focus2{1} = [0,0,0.6]; % Cartesian
    cancel2{1} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{1} = 1.2;

    focus2{2} = [0,0,0.6]; % Cartesian
    cancel2{2} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{2} = 1.2;

    focus2{3} = [0,0,0.6]; % Cartesian
    cancel2{3} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{3} = 1.2;

    focus2{4} = [0,0,0.6]; % Cartesian
    cancel2{4} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{4} = 1.2;

    focus2{5} = [0,0,0.6]; % Cartesian
    cancel2{5} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{5} = 1.2;

    focus2{6} = [0,0,0.6]; % Cartesian
    cancel2{6} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{6} = 1.2;

    focus2{7} = [0,0,0.6]; % Cartesian
    cancel2{7} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{7} = 1.2;

    focus2{8} = [0,0,0.6]; % Cartesian
    cancel2{8} = [0,0,0.4;0,0,0.76]; % Cartesian
    max_vol2{8} = 1.2;
    
    focus3{1} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{1} = [0,0,0.6]; % Cartesian
    max_vol3{1} = 1.7;

    focus3{2} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{2} = [0,0,0.6;]; % Cartesian
    max_vol3{2} = 1.7;

    focus3{3} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{3} = [0,0,0.6]; % Cartesian
    max_vol3{3} = 1.7;

    focus3{4} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{4} = [0,0,0.6]; % Cartesian
    max_vol3{4} = 1.7;

    focus3{5} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{5} = [0,0,0.6]; % Cartesian
    max_vol3{5} = 1.7;

    focus3{6} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{6} = [0,0,0.6]; % Cartesian
    max_vol3{6} = 1.7;

    focus3{7} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{7} = [0,0,0.6]; % Cartesian
    max_vol3{7} = 1.7;

    focus3{8} = [0,0,0.48;0,0,0.76]; % Cartesian
    cancel3{8} = [0,0,0.6]; % Cartesian
    max_vol3{8} = 1.7;

    patch = cell(n_patches, 1);
    cancel_mag1=num2cell(ones(n_patches,1));
    cancel_mag2=num2cell(ones(n_patches,1));
    cancel_mag3=num2cell(ones(n_patches,1));
    for j = 1:n_patches
	fprintf('Patch no. %d of %d\n', j, n_patches)
        phi_wrapped = wrapToPi(phi_patch{j});
        if phi_wrapped(1) > phi_wrapped(2)
            idx_phi = phi < phi_wrapped(2) | phi > phi_wrapped(1);
        else
            idx_phi = phi > phi_wrapped(1) & phi < phi_wrapped(2);
        end
        cnt = 1;

        while max(cancel_mag1{j})>0.95 || max(cancel_mag2{j})>0.95 || max(cancel_mag3{j})>0.9    % If cancellation is bad, randomly choose electrode locations again
        cnt = cnt+1;
        if cnt>5000
            error('Problem infeasible!')
        end
        patch{j} = surface(theta > theta_patch{j}(1) & theta < theta_patch{j}(2) & idx_phi,:);
        if size(patch{j}, 1) > num_pairs{j}
            rand_idx = randi(length(patch{j}), num_pairs{j}, 1);
            patch{j} = patch{j}(rand_idx, :);
        else
            num_pairs{j} = size(patch{j}, 1);
            warning('too many electrodes, using all possible ones.')
        end
        factors = zeros(num_nodes, num_pairs{j});
        for i = 1: num_pairs{j}
            [new_coord, mirror_x, mirror_y, swap_xy] = ConvertCoord(patch{j}(i,1),patch{j}(i,2),patch{j}(i,3), num_x, num_y);
            tmp = reshape(FactorTable(num2str(new_coord)), num_x, num_y, num_z);
            if swap_xy
                tmp = permute(tmp, [2,1,3]);
            end
            if mirror_x
                tmp = flip(tmp, 1);
            end
            if mirror_y
                tmp = flip(tmp, 2);
            end

            factors(:,i) = tmp(:);
        end
        [ElectrodeAmplitude1(:,j), currentmag1(:,j), cancel_mag1{j}, flag1] = spatial_interfere(focus1{j}, cancel1{j}, factors, max_vol1{j}, num_x, num_y, num_z);
        [ElectrodeAmplitude2(:,j), currentmag2(:,j), cancel_mag2{j}, flag2] = spatial_interfere(focus2{j}, cancel2{j}, factors, max_vol2{j}, num_x, num_y, num_z);
        [ElectrodeAmplitude3(:,j), currentmag3(:,j), cancel_mag3{j}, flag3] = spatial_interfere(focus3{j}, cancel3{j}, factors, max_vol3{j}, num_x, num_y, num_z);
        if flag1 || flag2 || flag3
            error('Problem infeasible!')
        end
        end
    %     disp(cancel_mag1{j})
    end
    fprintf('Finished optimization, starting HH simulation... take %d\n', take);
   
    %% HH Model
    CenterAmp1 = 107;
    CenterAmp2 = 107;
    CenterAmp3 = 110;
    f_center = 2000;      % center frequency
    delta_f = 30;
    if n_patches == 2
        freq_shift = [0;15];
    else
        freq_shift = [0:2*delta_f/n_patches:delta_f/2,delta_f/2-2*delta_f/n_patches:-2*delta_f/n_patches:-delta_f/2,-delta_f/2+2*delta_f/n_patches:2*delta_f/n_patches:-2*delta_f/n_patches]';
    end
    freq = f_center + freq_shift;

    Fs =  25; % kHz
    twidth = 1000; % ms
    vstart = -70; % mV

    idx_brain = find(nodes(:,3)>-0.01 & nodes(:,1).^2+nodes(:,2).^2+nodes(:,3).^2 <=r_brain^2+eps);

    fire_single = nan(size(currentmag1,1), 1);
    fire_double = fire_single;
    fire_other = fire_single;
    fire2_single = nan(length(idx_brain), 1);
    fire2_double = fire2_single;
    fire2_other = fire2_single;
    currentmag2_double = currentmag3(idx_brain, :);
    currentmag2_single = currentmag1(idx_brain, :);
    currentmag2_other = currentmag2(idx_brain, :);
    [b,a]=butter(3,[100/(Fs*500), 1000/(Fs*500)]);

    fprintf('Single focus...\n');
    parfor idx = 1:size(currentmag2_single, 1)
        amp = CenterAmp1 * currentmag2_single(idx,:)';
        [T,S] = HodgkinHuxleyFast(vstart,twidth,freq,amp,Fs);
        if length(S) < twidth*Fs+1
               fire2_single(idx) = -1;
            continue;
        end
        Vfilt = filtfilt(b,a,S);
        if max(Vfilt(10000:20000))>30
            fire2_single(idx) = 1;
        else
            fire2_single(idx) = 0;
        end
    end
    
    fprintf('Single focus done! Another single focus...\n');
    parfor idx = 1:size(currentmag2_other,1)
        amp = CenterAmp2 * currentmag2_other(idx,:)';
        [T,S] = HodgkinHuxleyFast(vstart,twidth,freq,amp,Fs);
        if length(S) < twidth*Fs+1
               fire2_other(idx) = -1;
            continue;
        end
        Vfilt = filtfilt(b,a,S);
        if max(Vfilt(10000:20000))>30
            fire2_other(idx) = 1;
        else
            fire2_other(idx) = 0;
        end
    end
    
    fprintf('Two foci...\n');
    parfor idx = 1:size(currentmag2_double,1)
        amp = CenterAmp3 * currentmag2_double(idx,:)';
        [T,S] = HodgkinHuxleyFast(vstart,twidth,freq,amp,Fs);
        if length(S) < twidth*Fs+1
               fire2_double(idx) = -1;
            continue;
        end
        Vfilt = filtfilt(b,a,S);
        if max(Vfilt(10000:20000))>30
            fire2_double(idx) = 1;
        else
            fire2_double(idx) = 0;
        end
    end
    
    

    fire_single(idx_brain) = fire2_single;
    fire_single = reshape(fire_single, num_x, num_y, num_z);

    fire_double(idx_brain) = fire2_double;
    fire_double = reshape(fire_double, num_x, num_y, num_z);
    
    fire_other(idx_brain) = fire2_other;
    fire_other = reshape(fire_other, num_x, num_y, num_z);
    
    toc;
    clear FactorTable
%     figure,slice(xgrid,ygrid,zgrid,permute(fire_single, [2,1,3]),0,0,[0.4])  % Visualize, change the last three numbers to the plane of interest
%     figure,slice(xgrid,ygrid,zgrid,permute(fire_double, [2,1,3]),0,0,[0.32])  % Visualize, change the last three numbers to the plane of interest
    save(['Results_corrected/STIMULUS_all80_',num2str(take)])
end
