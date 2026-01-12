%% TPMS HEAT SINK THERMAL MODEL (Correct + Extensive Comments + Fan h Estimate)
% Ambient: 35C everywhere (enforced)
% 2 CPUs feed the SAME copper (shared downstream path; corrected worst-junction constraint)
%
% What you get:
%   - Conduction resistances (chip/TIM/copper/TIM/base)
%   - h_min needed for each sink (Gyroid vs Schwarz-D), 1CPU and 2CPU
%   - Fan-based estimated h from 30 m^3/hr through your face areas (33x34mm, 58x34mm)
%   - Optional back-calculated "effective h" from Fusion Tmax values

clear; clc;
deg = char(176);

%% ========================== GIVEN MATERIAL PROPERTIES ==================
k_chip   = 150;   % W/m-K
k_TIM    = 10;    % W/m-K
k_copper = 400;   % W/m-K
k_HS     = 177;   % W/m-K

%% ========================== GIVEN GEOMETRY =============================
% Chip (die)
L_chip = 0.5e-3;                    % m
A_chip = (10e-3) * (10e-3);          % m^2

% TIM (both interfaces)
L_TIM_chip_to_copper = 0.1e-3;       % m
L_TIM_copper_to_HS   = 0.1e-3;       % m

% Copper spreader thickness
L_copper = 3e-3;                      % m

% HS base thickness
L_HS = 7e-3;                          % m

%% ========================== CONVECTION AREAS (Fusion) ==================
% Units: m^2
A_conv_gyroid_1CPU   = 15452.443e-6;
A_conv_gyroid_2CPU   = 29228.806e-6;

A_conv_schwarz_1CPU  = 19071.715e-6;
A_conv_schwarz_2CPU  = 32299.873e-6;

%% ========================== FOOTPRINT AREAS (Conduction) ===============
% 1 CPU footprint UPDATED per your input
A_copper_1CPU = (35e-3) * (36e-3);   % m^2  (35x36 mm)
A_HS_1CPU     = A_copper_1CPU;

% 2 CPU footprint (shared copper/base)
A_copper_2CPU = (60e-3) * (36e-3);   % m^2  (60x36 mm)
A_HS_2CPU     = A_copper_2CPU;

%% ========================== THERMAL LIMITS =============================
T_amb   = 35;       % degC
T_j_max = 85;       % degC
DeltaT_allow = T_j_max - T_amb; % degC

%% ========================== POWERS =====================================
P_1CPU = 20;        % W

P_CPU1 = 20;        % W
P_CPU2 = 15;        % W
P_tot  = P_CPU1 + P_CPU2;     % W
P_worst = max(P_CPU1, P_CPU2);

fprintf('========================================\n');
fprintf('INPUT SUMMARY\n');
fprintf('========================================\n');
fprintf('T_amb = %.1f %sC, T_j_max = %.1f %sC, DeltaT = %.1f %sC\n', T_amb, deg, T_j_max, deg, DeltaT_allow, deg);
fprintf('P_1CPU = %.1f W | 2CPU: P_tot = %.1f W (worst = %.1f W)\n', P_1CPU, P_tot, P_worst);

%% =======================================================================
%% 1 CPU: Conduction path + h_min
%% =======================================================================
% Network:
%   R_path_1 = R_chip + R_TIM(chip->copper)            [uses A_chip]
%   R_common_1 = R_copper + R_TIM(copper->HS) + R_HS   [uses footprint area]
%   R_conv = 1/(h*A_conv)
%
% Junction constraint:
%   Tj = Tamb + P*(R_noConv_1 + 1/(h*A_conv)) <= Tj_max
% => h_min = 1 / (A_conv*(DeltaT/P - R_noConv_1))

fprintf('\n========================================\n');
fprintf('1 CPU CALCULATIONS\n');
fprintf('========================================\n');

R_chip_1    = r_cond(L_chip, k_chip, A_chip);
R_TIM_CPU_1 = r_cond(L_TIM_chip_to_copper, k_TIM, A_chip);
R_path_1    = R_chip_1 + R_TIM_CPU_1;

R_copper_1  = r_cond(L_copper, k_copper, A_copper_1CPU);
R_TIM_HS_1  = r_cond(L_TIM_copper_to_HS, k_TIM, A_HS_1CPU);
R_HS_1      = r_cond(L_HS, k_HS, A_HS_1CPU);

R_common_1  = R_copper_1 + R_TIM_HS_1 + R_HS_1;
R_noConv_1  = R_path_1 + R_common_1;

fprintf('R_noConv_1 = %.6f %sC/W (chip+TIM+copper+TIM+base)\n', R_noConv_1, deg);

% Feasibility check
if (DeltaT_allow/P_1CPU) <= R_noConv_1
    error('1 CPU infeasible: conduction-only path already exceeds thermal budget.');
end

h_min_schwarz_1 = h_min_single_source(A_conv_schwarz_1CPU, P_1CPU, DeltaT_allow, R_noConv_1);
h_min_gyroid_1  = h_min_single_source(A_conv_gyroid_1CPU,  P_1CPU, DeltaT_allow, R_noConv_1);

fprintf('h_min (1CPU) Schwarz-D = %.3f W/m^2K\n', h_min_schwarz_1);
fprintf('h_min (1CPU) Gyroid    = %.3f W/m^2K\n', h_min_gyroid_1);

%% =======================================================================
%% 2 CPU: Correct worst-junction model for shared copper + base + convection
%% =======================================================================
% Topology:
%   Each CPU has its own R_path_2 into the copper (same as 1CPU path on chip area)
%   Then both feed shared downstream resistances R_common_2 and convection.
%
% Correct worst-junction constraint:
%   Tj_worst = Tamb + P_worst*R_path_2 + P_tot*(R_common_2 + 1/(h*A_conv)) <= Tj_max
%
% Solve:
%   1/(h*A_conv) <= (DeltaT - P_worst*R_path_2)/P_tot - R_common_2
%   h_min = 1/(A_conv * RHS)

fprintf('\n========================================\n');
fprintf('2 CPU CALCULATIONS (shared copper path)\n');
fprintf('========================================\n');

R_chip_2    = r_cond(L_chip, k_chip, A_chip);
R_TIM_CPU_2 = r_cond(L_TIM_chip_to_copper, k_TIM, A_chip);
R_path_2    = R_chip_2 + R_TIM_CPU_2;

R_copper_2  = r_cond(L_copper, k_copper, A_copper_2CPU);
R_TIM_HS_2  = r_cond(L_TIM_copper_to_HS, k_TIM, A_HS_2CPU);
R_HS_2      = r_cond(L_HS, k_HS, A_HS_2CPU);

R_common_2  = R_copper_2 + R_TIM_HS_2 + R_HS_2;

rhs_base = (DeltaT_allow - P_worst*R_path_2)/P_tot - R_common_2;
if rhs_base <= 0
    error('2 CPU infeasible: even infinite convection cannot meet junction limit.');
end

h_min_schwarz_2 = h_min_two_source_shared_path(A_conv_schwarz_2CPU, DeltaT_allow, P_worst, P_tot, R_path_2, R_common_2);
h_min_gyroid_2  = h_min_two_source_shared_path(A_conv_gyroid_2CPU,  DeltaT_allow, P_worst, P_tot, R_path_2, R_common_2);

fprintf('h_min (2CPU) Schwarz-D = %.3f W/m^2K\n', h_min_schwarz_2);
fprintf('h_min (2CPU) Gyroid    = %.3f W/m^2K\n', h_min_gyroid_2);

%% =======================================================================
%% FAN-BASED h ESTIMATE (from 30 m^3/hr through given face areas)
%% =======================================================================
% IMPORTANT:
% Fan "30 m^3/hr" is usually a free-air rating. Through a restrictive lattice,
% actual flow can drop significantly due to pressure drop.
%
% Therefore we include:
%   flow_factor (0..1): a derating factor applied to the fan nominal flow.
%   epsilon_open (0..1): fraction of face area that is actually open flow area.
%   Dh_model: a chosen hydraulic diameter for the flow passages.
%
% Why Dh_model matters:
%   If you use the full face dimensions as a duct, Dh is large and h comes out low.
%   In a TPMS lattice, the true channel Dh is much smaller (mm-scale), which increases h.

fprintf('\n========================================\n');
fprintf('FAN-BASED h ESTIMATE (engineering approximation)\n');
fprintf('========================================\n');

Q_fan_nom = 30/3600;     % m^3/s (30 m^3/hr)
flow_factor = 1.00;      % set <1 if you want to derate flow due to pressure drop (e.g., 0.5)
Q_eff = flow_factor * Q_fan_nom;

% Face areas (you gave these)
a1 = 33e-3; b1 = 34e-3;  % 1CPU face dimensions (m)
a2 = 58e-3; b2 = 34e-3;  % 2CPU face dimensions (m)

A_face_1 = a1*b1;
A_face_2 = a2*b2;

% Superficial velocities (based on FULL face area)
Vsup_1 = Q_eff / A_face_1;
Vsup_2 = Q_eff / A_face_2;

fprintf('Q_eff = %.6f m^3/s (flow_factor = %.2f)\n', Q_eff, flow_factor);
fprintf('1CPU face = %.0f x %.0f mm -> V_superficial = %.2f m/s\n', a1*1e3, b1*1e3, Vsup_1);
fprintf('2CPU face = %.0f x %.0f mm -> V_superficial = %.2f m/s\n', a2*1e3, b2*1e3, Vsup_2);

% Air properties around ~35C (approx constants; good enough for this level)
air.k  = 0.0263;     % W/m-K
air.Pr = 0.71;       % -
air.nu = 1.6e-5;     % m^2/s (kinematic viscosity)

% Open-area / porosity factor:
% If you do not know it, keep epsilon_open = 1.0 for a conservative (lower velocity) estimate.
% If you suspect ~60% open area, use epsilon_open = 0.6 (increases velocity, increases Re, increases h).
epsilon_open = 1.0;

% Hydraulic diameter model:
% OPTION A (very conservative): use face-as-duct hydraulic diameter
Dh_face_1 = rect_Dh(a1, b1);
Dh_face_2 = rect_Dh(a2, b2);

% OPTION B (more realistic for lattice): set a mm-scale Dh that represents typical flow channels
% If your cell size is ~4 mm, an effective Dh of ~3 mm is a reasonable starting point.
Dh_channel_guess = 3e-3;  % m  <-- change if you have better channel size

% Compute h using both models to give you a range.
[h1_face, Re1_face, Nu1_face, meth1] = h_from_flow(Q_eff, a1, b1, epsilon_open, Dh_face_1, air);
[h2_face, Re2_face, Nu2_face, meth2] = h_from_flow(Q_eff, a2, b2, epsilon_open, Dh_face_2, air);

[h1_chan, Re1_chan, Nu1_chan, meth1c] = h_from_flow(Q_eff, a1, b1, epsilon_open, Dh_channel_guess, air);
[h2_chan, Re2_chan, Nu2_chan, meth2c] = h_from_flow(Q_eff, a2, b2, epsilon_open, Dh_channel_guess, air);

fprintf('\n--- Using face-as-duct Dh (very conservative) ---\n');
fprintf('1CPU: Dh=%.2f mm, Re=%.0f, Nu=%.1f (%s) -> h≈%.2f W/m^2K\n', Dh_face_1*1e3, Re1_face, Nu1_face, meth1, h1_face);
fprintf('2CPU: Dh=%.2f mm, Re=%.0f, Nu=%.1f (%s) -> h≈%.2f W/m^2K\n', Dh_face_2*1e3, Re2_face, Nu2_face, meth2, h2_face);

fprintf('\n--- Using channel Dh guess (more realistic for TPMS) ---\n');
fprintf('1CPU: Dh=%.2f mm, Re=%.0f, Nu=%.1f (%s) -> h≈%.2f W/m^2K\n', Dh_channel_guess*1e3, Re1_chan, Nu1_chan, meth1c, h1_chan);
fprintf('2CPU: Dh=%.2f mm, Re=%.0f, Nu=%.1f (%s) -> h≈%.2f W/m^2K\n', Dh_channel_guess*1e3, Re2_chan, Nu2_chan, meth2c, h2_chan);

% Predict temperatures using estimated h (pick the model you want to use)
% Here I use the channel-based estimate as the default for lattice flow.
h_use_1 = h1_chan;
h_use_2 = h2_chan;

Tj_schwarz_1_pred = T_amb + P_1CPU * (R_noConv_1 + 1/(h_use_1 * A_conv_schwarz_1CPU));
Tj_gyroid_1_pred  = T_amb + P_1CPU * (R_noConv_1 + 1/(h_use_1 * A_conv_gyroid_1CPU));

Tj_schwarz_2_pred = T_amb + P_worst*R_path_2 + P_tot*(R_common_2 + 1/(h_use_2 * A_conv_schwarz_2CPU));
Tj_gyroid_2_pred  = T_amb + P_worst*R_path_2 + P_tot*(R_common_2 + 1/(h_use_2 * A_conv_gyroid_2CPU));

fprintf('\n--- Predicted junction temps using h_use (channel Dh guess) ---\n');
fprintf('1CPU Schwarz-D: Tj_pred = %.2f %sC (limit %.1f)\n', Tj_schwarz_1_pred, deg, T_j_max);
fprintf('1CPU Gyroid:    Tj_pred = %.2f %sC (limit %.1f)\n', Tj_gyroid_1_pred,  deg, T_j_max);
fprintf('2CPU Schwarz-D: Tj_pred = %.2f %sC (limit %.1f)\n', Tj_schwarz_2_pred, deg, T_j_max);
fprintf('2CPU Gyroid:    Tj_pred = %.2f %sC (limit %.1f)\n', Tj_gyroid_2_pred,  deg, T_j_max);

%% =======================================================================
%% SUMMARY
%% =======================================================================
fprintf('\n========================================\n');
fprintf('SUMMARY - h_min (required)\n');
fprintf('========================================\n');
fprintf('1 CPU: Schwarz-D %.3f | Gyroid %.3f  (W/m^2K)\n', h_min_schwarz_1, h_min_gyroid_1);
fprintf('2 CPU: Schwarz-D %.3f | Gyroid %.3f  (W/m^2K)\n', h_min_schwarz_2, h_min_gyroid_2);

fprintf('\nDone.\n');
%% =======================================================================
%% DASHBOARD PLOTS (Clean, single window)  --- UPDATED FUSION Tmax LIST ---
%% =======================================================================

cfgNames  = {'1 CPU','2 CPU'};
sinkNames = {'Gyroid','Schwarz-D'};

% Convection areas in mm^2 (so no scientific notation needed)
Aconv_mm2 = [A_conv_gyroid_1CPU,  A_conv_schwarz_1CPU; ...
             A_conv_gyroid_2CPU,  A_conv_schwarz_2CPU] * 1e6;

% h_min matrix (W/m^2K): rows=config, cols=sink
hmin = [h_min_gyroid_1,  h_min_schwarz_1; ...
        h_min_gyroid_2,  h_min_schwarz_2];

% ---- DeltaT breakdown at h_min (stacked) ----
% 1CPU: DeltaT = P*(R_noConv_1 + 1/(h*A))
dT1_cond    = P_1CPU * R_noConv_1;
dT1_conv_gy = P_1CPU * (1/(h_min_gyroid_1  * A_conv_gyroid_1CPU));
dT1_conv_sc = P_1CPU * (1/(h_min_schwarz_1 * A_conv_schwarz_1CPU));

% 2CPU worst junction: DeltaT = P_worst*R_path_2 + P_tot*(R_common_2 + 1/(h*A))
dT2_local       = P_worst * R_path_2;
dT2_shared_cond = P_tot   * R_common_2;
dT2_conv_gy     = P_tot * (1/(h_min_gyroid_2  * A_conv_gyroid_2CPU));
dT2_conv_sc     = P_tot * (1/(h_min_schwarz_2 * A_conv_schwarz_2CPU));

caseNames = {'1CPU Gy','1CPU Sch','2CPU Gy','2CPU Sch'};

% Columns = [LocalPath, SharedConduction, SharedConvection]
% For 1CPU, we place conduction into SharedConduction and LocalPath=0 for consistent legend.
dT_stack = [ ...
    0,         dT1_cond,        dT1_conv_gy; ...
    0,         dT1_cond,        dT1_conv_sc; ...
    dT2_local, dT2_shared_cond, dT2_conv_gy; ...
    dT2_local, dT2_shared_cond, dT2_conv_sc];

% ---- Fusion Tmax (max anywhere in solid) ----
% NOTE: You said Gyroid 1CPU no-fan is not done yet -> set as NaN (blank bar).
Tmax_fusion = [ ...
    47.27;     % Gyroid 1CPU with fan
    165.614;       % Gyroid 1CPU no fan (PENDING)
    79.316;    % Gyroid 2CPU with fan
    191.817;   % Gyroid 2CPU no fan
    46.629;    % Schwarz-D 1CPU with fan
    161.853;   % Schwarz-D 1CPU no fan
    54.811;    % Schwarz-D 2CPU with fan
    199.002];  % Schwarz-D 2CPU no fan

fusionShortLabels = { ...
    'Gy 1 Fan','Gy 1 NoFan','Gy 2 Fan','Gy 2 NoFan', ...
    'Sch 1 Fan','Sch 1 NoFan','Sch 2 Fan','Sch 2 NoFan'};

% ================== Figure ==================
fig = figure('Name','TPMS Summary Dashboard','Position',[80 80 1400 800]);
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(t,'TPMS Heat Sink Summary');

% Consistent axes styling
setAxesNice = @(ax) set(ax,'FontName','Calibri','FontSize',11,'LineWidth',1,'Box','on');

%% (1) A_conv comparison
ax1 = nexttile;
hb1 = bar(categorical(cfgNames), Aconv_mm2, 'grouped');
ylabel('A_{conv} (mm^2)');
title('Effective Convection Area (Fusion)');
grid on;
setAxesNice(ax1);

% Force normal tick display (no ×10^n)
ax1.YAxis.Exponent = 0;

legend(hb1, sinkNames, 'Location','northoutside','Orientation','horizontal');

% Value labels (safe placement using XEndPoints if available)
for j = 1:numel(hb1)
    if isprop(hb1(j),'XEndPoints')
        x = hb1(j).XEndPoints; y = hb1(j).YEndPoints;
        text(x, y, compose('%.0f',y), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize',9);
    end
end

%% (2) h_min comparison
ax2 = nexttile;
hb2 = bar(categorical(cfgNames), hmin, 'grouped');
ylabel('h_{min} (W/m^2K)');
title('Minimum Required Convection Coefficient');
grid on;
setAxesNice(ax2);

legend(hb2, sinkNames, 'Location','northoutside','Orientation','horizontal');

for j = 1:numel(hb2)
    if isprop(hb2(j),'XEndPoints')
        x = hb2(j).XEndPoints; y = hb2(j).YEndPoints;
        text(x, y, compose('%.1f',y), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize',9);
    end
end

%% (3) DeltaT budget (stacked)
ax3 = nexttile;
hb3 = bar(categorical(caseNames), dT_stack, 'stacked');
ylabel(['\DeltaT contributions (' char(176) 'C)']);
title('Temperature-Rise Budget at h_{min} (sums to 50°C)');
grid on;
setAxesNice(ax3);

legend(hb3, {'Local Path','Shared Conduction','Shared Convection'}, ...
    'Location','northoutside','Orientation','horizontal');

% Add the 50°C limit line but hide it from the legend
hold on;
if exist('yline','file') == 2
    yl = yline(DeltaT_allow,'--','LineWidth',1.5);
    yl.HandleVisibility = 'off';
else
    xl = xlim;
    pl = plot(xl,[DeltaT_allow DeltaT_allow],'--','LineWidth',1.5);
    set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold off;

ylim([0, max(DeltaT_allow*1.15, max(sum(dT_stack,2))*1.15)]);

%% (4) Fusion Tmax vs limit (UPDATED + supports NaN blank bar)
ax4 = nexttile;

catLabels = categorical(fusionShortLabels);
hbar4 = bar(catLabels, Tmax_fusion);
ylabel(['Fusion T_{max} (' char(176) 'C)']);
title('Fusion T_{max} (Max anywhere in solid) vs Limit');
grid on;
setAxesNice(ax4);

% Temperature limit line (excluded from legend)
hold on;
if exist('yline','file') == 2
    yl2 = yline(T_j_max,'--','LineWidth',1.8);
    yl2.HandleVisibility = 'off';
else
    xl = xlim;
    pl2 = plot(xl,[T_j_max T_j_max],'--','LineWidth',1.8);
    set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
hold off;

% Make sure y-limits include the largest non-NaN value
Tmax_max = max(Tmax_fusion(~isnan(Tmax_fusion)));
ylim([0, max(T_j_max*1.15, Tmax_max*1.05)]);

% PASS/FAIL labels (skip NaN; show PENDING)
for i = 1:numel(Tmax_fusion)
    if isnan(Tmax_fusion(i))
        text(i, 5, "PENDING", 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',9);
        continue;
    end
    status = "PASS";
    if Tmax_fusion(i) > T_j_max, status = "FAIL"; end

    text(i, Tmax_fusion(i), " " + status, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',9);
end

% Note inside the tile (won't overlap x-labels)
text(ax4, 0.02, 0.02, ...
    'Note: Fusion T_{max} is max anywhere in solid (not necessarily chip/junction).', ...
    'Units','normalized','FontSize',9,'VerticalAlignment','bottom');

%% ============================= FUNCTIONS ===============================

function R = r_cond(L, k, A)
% r_cond  1D conduction thermal resistance of a uniform slab.
%   R = L/(k*A)
% Inputs:
%   L [m], k [W/m-K], A [m^2]
% Output:
%   R [K/W] (numerically same as °C/W)
    R = L / (k * A);
end

function hmin = h_min_single_source(A_conv, P, DeltaT_allow, R_noConv)
% h_min_single_source  Minimum convection coefficient for single heat source.
% Model:
%   Tj = Tamb + P*(R_noConv + 1/(h*A_conv))
% Constraint:
%   Tj <= Tj_max  ->  1/(h*A_conv) <= DeltaT/P - R_noConv
% => h_min = 1 / (A_conv*(DeltaT/P - R_noConv))
    hmin = 1 / (A_conv * (DeltaT_allow/P - R_noConv));
end

function hmin = h_min_two_source_shared_path(A_conv, DeltaT_allow, P_worst, P_tot, R_path, R_common)
% h_min_two_source_shared_path  Minimum h for 2 sources feeding shared path.
% Worst junction model:
%   Tj_worst = Tamb + P_worst*R_path + P_tot*(R_common + 1/(h*A_conv))
% Constraint:
%   Tj_worst <= Tj_max
% => 1/(h*A_conv) <= (DeltaT - P_worst*R_path)/P_tot - R_common
    rhs = (DeltaT_allow - P_worst*R_path)/P_tot - R_common;
    hmin = 1 / (A_conv * rhs);
end

function Dh = rect_Dh(a, b)
% rect_Dh  Hydraulic diameter for a rectangular duct of width a and height b:
%   Dh = 2ab/(a+b)
    Dh = 2*a*b/(a+b);
end

function [h, Re, Nu, method] = h_from_flow(Q, a, b, epsilon_open, Dh_model, air)
% h_from_flow  Estimate convection coefficient h from volumetric flow Q through a face.
%
% Steps:
%   1) Compute FACE area A_face = a*b
%   2) Convert to "interstitial" velocity using open-area factor epsilon_open:
%        V = Q / (epsilon_open*A_face)
%      If epsilon_open < 1, actual velocity in the open channels increases.
%   3) Compute Reynolds number:
%        Re = V*Dh/nu
%   4) Use a basic correlation to estimate Nusselt number Nu:
%        - Laminar (Re < 2300): Nu = 3.66 (fully developed, constant wall temp)
%        - Turbulent (Re >= 10000): Dittus-Boelter Nu = 0.023 Re^0.8 Pr^0.4
%        - Transitional: use Dittus-Boelter but label as "TRANSITIONAL"
%   5) Convert to h:
%        h = Nu*k/Dh
%
% IMPORTANT LIMITATION:
%   This treats the lattice as an "equivalent duct/wall" problem. In a TPMS porous
%   lattice, the true Dh and open area can differ significantly, and pressure drop
%   can reduce Q. This is an engineering estimate, not a CFD replacement.

    A_face = a*b;
    V = Q / (epsilon_open * A_face);

    Re = V * Dh_model / air.nu;

    if Re < 2300
        Nu = 3.66;
        method = 'LAMINAR Nu=3.66';
    elseif Re >= 10000
        Nu = 0.023 * (Re^0.8) * (air.Pr^0.4);
        method = 'TURB Dittus-Boelter';
    else
        Nu = 0.023 * (Re^0.8) * (air.Pr^0.4);
        method = 'TRANSITIONAL (DB used)';
    end

    h = Nu * air.k / Dh_model;
end
