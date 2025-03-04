%% 
% Initial Parameters

L = 1200;     % Length of bridge
n = 1200;     % Discretize into 1 mm seg.
P = 452;      % Total weight of train [N]
x = linspace(0, L, n+1);     % x-axis
tensile_s = 30; % Mpa
compressive_s = 6; % Mpa
material_shear_strength = 4; % MPa
E = 4000; % MPa
mu = 0.2; % Poisson's ratio
glue_shear_strength = 2; % MPa 
%% 
% Train Details

train_length = 960; % Train length [mm]
x_train = [52 228 392 568 732 908]; % Train Axle Locations


d_axle1 =  x_train(1);
e_train_length = x_train(6) - d_axle1;
e_x_train = x_train - d_axle1; % Assuming loading begins when first axle reaches bridge
n_train = n - e_train_length; % Train can start from -52mm to 240mm
%% 
% Load Case 0

%{
P_axle = P / 6; % Load per axle [N]
P_train = [1 1 1 1 1 1] * P/6;  
%}
%% 
% Load Case 1


% Base Case: Freight cars have identical weights
freight_weight_base = P / (2 + 1.35); % Each freight car weight
locomotive_weight_base = 1.35 * freight_weight_base;

P_train = [freight_weight_base/2, freight_weight_base/2, ...   
freight_weight_base/2, freight_weight_base/2, ...
locomotive_weight_base/2, locomotive_weight_base/2,];
%}
%% 
% Load Case 2

%{ 
%Weight ratios
heavier_to_lighter_ratio = 1.10;
locomotive_to_heaviest_ratio = 1.35;
min_weight = P / (1 + heavier_to_lighter_ratio + locomotive_to_heaviest_ratio * heavier_to_lighter_ratio);

lighter_freight_weight = min_weight;
heavier_freight_weight = heavier_to_lighter_ratio * min_weight;
locomotive_weight = locomotive_to_heaviest_ratio * heavier_freight_weight;

% Distribute weights across axles
P_train = [heavier_freight_weight/2, heavier_freight_weight/2, ...
           lighter_freight_weight/2, lighter_freight_weight/2, ... % Lighter freight car axles
           locomotive_weight/2, locomotive_weight/2]; % Heavier freight car axles
%}
%% 
% SFD and BMD Conditions

SFDi = zeros(n_train+1, n+1);      % 1 SFD for each train loc.
BMDi = zeros(n_train+1, n+1);      % 1 BMD for each train loc.
SFDi(1, :) = x;
BMDi(1, :) = x;
%% 
% Find the Largest Possible Shear Force and Bending Moment for Different Locations

for iterations = 0:n_train

    % Start location of the train
    train_start = iterations;
    train_axles = train_start + e_x_train; % Axle positions on bridge

    % Initialize the load distribution matrix
    bridge_load = zeros(2, n+1); % 2 rows: 1st row for positions, 2nd row for loads

    % Fill the first row with the bridge positions (0mm, 1mm, 2mm, ...1199mm, 1200mm)
    bridge_load(1, :) = 0:n; 

    % Assign loads to the corresponding positions of the axles
    for axle = 1:length(train_axles)
        bridge_load(2, train_axles(axle) + 1) = P_train(axle); % Add axle load to corresponding position
    end

    % Reaction forces
    % Sum of moments at A = 0
    Reaction_B = sum(bridge_load(2, :) .* x) / L;
    % Sum of vertical Forces = 0
    Reaction_A = sum(P - Reaction_B);
    
    % Start with the reaction at A
    SFDi(iterations+1, 1) = Reaction_A;

    % Calculate shear force at each point
    for j = 2:n+1
        % Shear = Previous shear - Load at this point
        SFDi(iterations+1, j) = SFDi(iterations+1, j-1) - bridge_load(2, j); 
    end

    for j = 2:n+1
        % BMD at position j = Previous BMD + Shear Force at j-1 * Segment Length
        BMDi(iterations+1, j) = BMDi(iterations+1, j-1) + (SFDi(iterations+1, j-1)/1000);
    end
    BMDi(:, 1201) = 0;
end
%% 
% Get The SFD and BMD envelopes

SFD = max(abs((SFDi(2:end, :))));
BMD = max((BMDi(2:end, :))); 

BMD_subset = BMDi(2:end, :);  % Exclude the first row
SFD_subset = SFDi(2:end, :);  % Exclude the first row

% Find the indices of the maximum magnitude values
[~, idxBMD] = max(abs(BMD_subset), [], 1);
[~, idxSFD] = max(abs(SFD_subset), [], 1);

% Use the indices to retrieve the signed values
BMD_envelope = BMD_subset(sub2ind(size(BMD_subset), idxBMD, 1:size(BMD_subset, 2)));
SFD_envelope = SFD_subset(sub2ind(size(SFD_subset), idxSFD, 1:size(SFD_subset, 2)));

BMDmax = max(BMD)
SFDmax = max(SFD)
[max_BMD_row, ~] = find(abs(BMDi) == BMDmax, 1, 'first'); % Find row of max BMD
[max_SFD_row, ~] = find(abs(SFDi) == SFDmax, 1, 'first'); % Find row of max BMD
first_axle_position_for_BMDmax = max_BMD_row - 1
first_axle_position_for_SFDmax = max_SFD_row - 1
%% 
% Bridge Parameters

% = xc, tfw, tft, gcw, gch, rwh, lwh, wt, bfw, bft
% Design 0 Case
extra_h_tf = 1.27;
extra_h_webs = 0;
extra_b_gluetabs = 20;
x_c = (0:1200)'; 
tfb = 100 * ones(size(x_c)); 
tfh = (1.27 + extra_h_tf) * ones(size(x_c)); 
lgcb = (6.27 + extra_b_gluetabs) * ones(size(x_c));
lgch = 1.27 * ones(size(x_c));
rgcb = (6.27 + extra_b_gluetabs) * ones(size(x_c));
rgch = 1.27 * ones(size(x_c));
rwh = (114.92 + extra_h_webs) * ones(size(x_c));
lwh = (114.92 + extra_h_webs) * ones(size(x_c));
rwb = 1.27 * ones(size(x_c));
lwb = 1.27 * ones(size(x_c));
bfb = 80 * ones(size(x_c));
bfh = 1.27 * ones(size(x_c));

param = [x_c, tfb, tfh, lgcb, rgcb, lgch, rgch, rwh, lwh, rwb, lwb, bfb, bfh]; % Combine into a single matrix


% tfb = interp1(param(:,1), param(:,2), x);
% tfh = interp1(param(:,1), param(:,3), x);
% x_c Location, x, of cross-section change
% tfb Top Flange Width
% tfh Top Flange Thickness
% lgcb Left Glue Connector Width
% rgcb Right Glue Connector Width
% lgch Left Glue Connector Height
% rgch Right Glue Connector Height
% rwh Right Web Height
% lwh Left Web Height
% rwb Right Web Thickness
% lwb Left Web Thickness
% btb Bottom Flange Width
% bfh Bottom Flange Thickness
%% 
% Create Shape with Polyin

bs = [bfb, lwb, rwb, lgcb, rgcb, tfb];
sorted_bs = sort(bs, 2, 'descend');
total_b = sorted_bs(:, 1);
cs_x = 1;
total_h = (bfh + lwh + lgch + tfh);

side = (total_b - sorted_bs(:, 2)) / 2;
% Top Flange
x_tf = [zeros(size(x_c)) tfb tfb zeros(size(x_c))];
y_tf = [total_h-tfh total_h-tfh total_h total_h];
csTop_flange = polyshape({x_tf(cs_x,:)}, {y_tf(cs_x,:)});
% Left and Right Webs
x_wL = [side (side + lwb) (side + lwb) side];
y_wL = [bfh bfh (bfh + lwh) (bfh + lwh)];
csLeft_web = polyshape({x_wL(cs_x,:)}, {y_wL(cs_x,:)});
x_wR = [(total_b - side - rwb) (total_b - side) (total_b - side) (total_b - side - rwb)];
y_wR = [bfh bfh (bfh + lwh) (bfh + lwh)];
csRight_web = polyshape({x_wR(cs_x,:)}, {y_wR(cs_x,:)});
% Left and Right Glue Connectors
x_gcL = [side (side + lgcb) (side + lgcb) side];
y_gcL = [(bfh + lwh) (bfh + lwh) (bfh + lwh + lgch) (bfh + lwh + lgch)];
csLeft_glue_connector = polyshape({x_gcL(cs_x,:)}, {y_gcL(cs_x,:)});
x_gcR = [(total_b - side - rgcb) (total_b - side) (total_b - side) (total_b - side - rgcb)];
y_gcR = [(bfh + lwh + rgch) (bfh + lwh + rgch) (bfh + lwh) (bfh + lwh)];
csRight_glue_connector = polyshape({x_gcR(cs_x,:)}, {y_gcR(cs_x,:)});
% Bottom Flange
x_bf = [side (total_b - side) (total_b - side) side];
y_bf = [zeros(size(x_c)) zeros(size(x_c)) bfh bfh];
csBottom_flange = polyshape({x_bf(cs_x,:)}, {y_bf(cs_x,:)});
% Example Cross-section
all_cs = [csTop_flange csLeft_web csRight_web csLeft_glue_connector csRight_glue_connector csBottom_flange];

% Iteratively union the polyshapes
beam_cs  = all_cs(1);
for i = 2:length(all_cs)
    beam_cs  = union(beam_cs, all_cs(i));
end
%% 
% Find Centroids of the Cross-Section and Individual Components

[x_centroid, y_centroid] = centroid(beam_cs)
figure;
plot(beam_cs)
hold on
plot(x_centroid, y_centroid, "r*")
hold off


[tempx, tempy] = centroid(csTop_flange);
tf_centroid = [tempx tempy];
[tempx, tempy] = centroid(csLeft_web);
lw_centroid = [tempx tempy];
[tempx, tempy] = centroid(csRight_web);
rw_centroid = [tempx tempy];
[tempx, tempy] = centroid(csLeft_glue_connector);
lgc_centroid = [tempx tempy];
[tempx, tempy] = centroid(csRight_glue_connector);
rgc_centroid = [tempx tempy];
[tempx, tempy] = centroid(csBottom_flange);
bf_centroid =[tempx tempy];

%{
Design 0 Dimensions:
rgcb = 6.27
lgcb = 6.27
rgch = 1.27
lgch = 1.27
rwh = 72.46
lwh = 72.46
rwb = 1.27
lwb = 1.27
bfb = 80
bfh = 1.27
%}
%% 
% Find the Area of the Components

A_tf = area(csTop_flange);
A_lw = area(csLeft_web);
A_rw = area(csRight_web);
A_lgc = area(csLeft_glue_connector);
A_rgc = area(csRight_glue_connector);
A_bf = area(csBottom_flange);
%% 
% Find the Distance of Each Component to the Centroid

% Centroid of Component - Cross-section Centroid
d_tf = tf_centroid(2) - y_centroid;
d_lw = lw_centroid(2) - y_centroid;
d_rw = rw_centroid(2) - y_centroid;
d_lgc = lgc_centroid(2) - y_centroid;
d_rgc = rgc_centroid(2) - y_centroid;
d_bf = bf_centroid(2) - y_centroid;
%% 
% Find the Moment of Area of each Component

% bh^3/12
I_tf = tfb(cs_x) * tfh(cs_x)^3 / 12;
I_lw = lwb(cs_x) * lwh(cs_x)^3 / 12;
I_rw = rwb(cs_x) * rwh(cs_x)^3 / 12;
I_lgc = lgcb(cs_x) * lgch(cs_x)^3 / 12;
I_rgc = rgcb(cs_x) * rgch(cs_x)^3 / 12;
I_bf = bfb(cs_x) * bfh(cs_x)^3 / 12;

%% 
% Find the Moment of Area of the Cross-Section

I_yy = ... 
    (I_tf + A_tf * (d_tf)^2) + ...
    (I_lw + A_lw * (d_lw)^2) + ...
    (I_rw + A_rw * (d_rw)^2) + ...
    (I_lgc + A_lgc * (d_lgc)^2) + ...
    (I_rgc + A_rgc * (d_rgc)^2) + ...
    (I_bf + A_bf * (d_bf)^2)
%% 
% Sectional Properties

ybar = y_centroid;
ybot = abs(y_centroid - bf_centroid(2));
ytop = abs(tf_centroid(2) - y_centroid);
% Q at centroidal axes
A_cent1 = (y_centroid - bfh(cs_x)) * lwb(cs_x);
A_cent2 = (y_centroid - bfh(cs_x)) * rwb(cs_x);
d_cent1 = (y_centroid - bfh(cs_x)) / 2;
d_cent2 = (y_centroid - bfh(cs_x)) / 2;
Qcent = A_bf * abs(d_bf) + A_cent1 * abs(d_cent1) + A_cent2 * abs(d_cent2)
bcent = lwb(cs_x) + rwb(cs_x);
% Q at glue location
glue_widthL = lgcb(cs_x) - lwb(cs_x); % Glue width on left glue connector
glue_widthR = rgcb(cs_x) - rwb(cs_x); % Glue width on right connector
Qglue = A_tf * d_tf
bglue = glue_widthL + glue_widthR;
%% 
% Calculate Applied Stress

% Use the Largest Bending Moment and Shear Force
S_top = BMDmax * 1000 * (total_h(cs_x) - y_centroid) / I_yy
S_bot = BMDmax * 1000 * y_centroid / I_yy
shear_cent = (SFDmax * Qcent) / (I_yy * bcent)
shear_glue = (SFDmax * Qglue) / (I_yy * bglue)
%% 
% Calculate Other Failure Modes
%% Material and Thin Plate Buckling Capacities

diaphram_distance = 115; % Distance Between Diaphrams
S_tens = S_bot;
S_comp = S_top;
shear_max = shear_cent;
shear_gmax = shear_glue;
% Buckling of compressive flange between the webs
S_buck1 = (4 * pi^2 * E)/(12 * (1 - mu^2)) * ((tfh(cs_x) / (bfb(cs_x) - lwb(cs_x)/2 - rwb(cs_x)/2))^2);
% Buckling of Flange Tips
S_buck2 = ((0.425 * pi^2 * E)/(12 * (1 - mu^2))) * ((tfh(cs_x) / (((tfb(cs_x) - bfb(cs_x))/2) + rwb(cs_x)/2)))^2;
% Buckling of webs due to variable flexural stresses
% This is for the left web but for simplicity assume webs have same thickness
S_buck3 = (6 * pi^2 * E)/(12 * (1 - mu^2)) * ((lwb(cs_x) / (total_h(cs_x) - y_centroid))^2);
% For Right web if right web hass different thickness:
%S_buck3R = (6 * pi^2 * E)/(12 * (1 - mu^2)) * ((rwb(cs_x) / (total_h(cs_x) - tfh(cs_x) - y_centroid))^2)
% Shear Buckling of the webs
shear_buck = ((5 * pi^2 * E)/(12*(1-mu^2))) * ((lwb(cs_x)/(total_h(cs_x) - tfh(cs_x) - bfh(cs_x)))^2 + ...
    (lwb(cs_x)/diaphram_distance)^2);
%% 
% Factor of Safeties

FOS_tens = tensile_s / S_tens
FOS_comp = compressive_s / S_comp
FOS_shear = material_shear_strength / shear_max
FOS_glue = glue_shear_strength / shear_gmax
FOS_buck1 = S_buck1 / S_comp
FOS_buck2 = S_buck2 / S_comp
FOS_buck3 = S_buck3 / S_comp
FOS_buckV = shear_buck / ((SFDmax * Qcent) / (I_yy * (rwb(cs_x) + lwb(cs_x))))
%% 
% Minimum Factor of Safety and Load

minFOS = min([FOS_tens, FOS_comp, FOS_shear, FOS_glue, FOS_buck1, FOS_buck2, FOS_buck3, FOS_buckV])
Pf = P * minFOS
%% 
% Calculating Failure Loads

Mf_tens = BMDmax * FOS_tens;
Mf_comp = BMDmax * FOS_comp;

Vf_shear = (material_shear_strength * I_yy * bcent)/Qcent;
Vf_glue = ((glue_shear_strength * I_yy * bglue)/Qglue);
Mf_buck1 = (S_buck1 * I_yy) / (total_h(cs_x) - y_centroid) / 1000
Mf_buck2 = (S_buck2 * I_yy) / (total_h(cs_x) - y_centroid) / 1000
Mf_buck3 = (S_buck3 * I_yy) / (total_h(cs_x) - y_centroid) / 1000
Vf_buckV = (shear_buck * I_yy * bcent) / Qcent;
%% 
% Plotting

% Shear Force Failures
figure;
hold off;
subplot(2,3,1);
hold on; grid on; grid minor;
plot(x, Vf_shear * ones(1, L + 1), 'r', 'DisplayName', 'Matboard Shear Failure');
plot(x, -Vf_shear * ones(1, L + 1), 'r', 'HandleVisibility', 'off');
plot(x, SFD_envelope, 'k', 'HandleVisibility', 'off');
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
legend;
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');

subplot(2,3,2);
title("Shear Force Diagrams vs Shear Force Capacities")
hold on; grid on; grid minor;
plot(x, Vf_glue * ones(1, L + 1), 'r', 'DisplayName', 'Glue Shear Failure')
plot(x, -Vf_glue * ones(1, L + 1), 'r', 'HandleVisibility', 'off')
plot(x, SFD_envelope, 'k', 'HandleVisibility', 'off')
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
legend();
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');


subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, Vf_buckV * ones(1, L + 1), 'r', 'DisplayName', 'Matboard Shear Failure')
plot(x, -Vf_buckV * ones(1, L + 1), 'r', 'HandleVisibility', 'off')
plot(x, SFD_envelope, 'k', 'HandleVisibility', 'off')
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
legend();
xlabel('Distance along bridge (mm)');
ylabel('Shear Force (N)');

% Bending Moment Failures

subplot(2,3,4);
hold on; grid on; grid minor;
plot(x, Mf_tens * ones(1, L + 1), 'r', 'DisplayName', 'Matboard Tension Failure')
plot(x, Mf_comp * ones(1, L + 1), 'b', 'DisplayName', "Matboard Compression Failure")
plot(x, BMD_envelope, 'k', 'HandleVisibility', 'off')
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
ylim([-50, max([BMDmax, Mf_tens, Mf_comp]) * 1.2]); 
set(gca, 'YDir', 'reverse');
legend('Location', 'southeast');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nm)');

subplot(2,3,5)
title("Bending Moment Diagrams vs Bending Moment Capacities")
hold on; grid on; grid minor;
plot(x, Mf_buck1 * ones(1, L + 1), 'r', 'DisplayName', "Matboard Buckling Failure, Top Flange, Mid")
plot(x, Mf_buck2 * ones(1, L + 1), 'b', 'DisplayName', 'Matboard Buckling Failure, Top Flange, Sides')
plot(x, BMD_envelope, 'k', 'HandleVisibility', 'off')
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
ylim([-50, max([BMDmax, Mf_buck1, Mf_buck2]) * 1.2]); 
set(gca, 'YDir', 'reverse');
legend('Location', 'southeast');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nm)');

subplot(2,3,6);
hold on; grid on; grid minor;
plot(x, Mf_buck3 * ones(1, L + 1), 'r', 'DisplayName', 'Matboard Buckling Failure, Webs')
plot(x, BMD_envelope, 'k', 'HandleVisibility', 'off')
plot([0, L], [0, 0], 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
ylim([-50, max([BMDmax, Mf_buck3]) * 1.2]); 
set(gca, 'YDir', 'reverse');
legend('Location', 'southeast');
xlabel('Distance along bridge (mm)');
ylabel('Bending Moment (Nm)');
hold off;

[~, maxBMD_index] = max(BMD_envelope(1, :));
maxBMD_location = maxBMD_index - 1
csA_total = (A_tf + A_rgc + A_lgc + A_lw + A_rw + A_bf)
Volume_total = csA_total * 1260