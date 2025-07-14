function [result, t] = run_hiMCES(args)
% run_hiMCES   Run hiMCES
%   Version Disease Models & Mechanisms (https://doi.org/10.1242/dmm.050365)
%
%   This method has named arguments with default values
%   NAME = DEFAULT  description
%
%		Simulation parameters
%    simTime = 800       Simulation time
%    MaxStep = 1e-3      Longest allowed step for ode
%    InitialStep = 2e-5  Initial step for ode
%
%   Output parameters
%    result = 'Vm'  Values stored into the result defined with comma separated list, in order
%
%   Stimulation settings
%    stimFlag = 0	 Stimulation OFF:0 / ON:1
%
%   Drug dosage
%    tDrugApplication = 1400  Time to start applying the drug
%    dose_uM = 0.3            Drug dosage (uM), allowed values are: 0.3 / 2 / 10
%    INaLRedMed = 1           OFF:0 / ON:1
%    IKsRedMed = 1            OFF:0 / ON:1
%    INaCaRedMed = 1          OFF:0 / ON:1
%    IfRedMed = 1             OFF:0 / ON:1
%
%   Cell state
%    ischemiaFlag = 0        no ischemia:0 / SEV1: 1 / SEV2: 2
%    HyperkalemiaFlag = 0    Hyperkalemia OFF:0 / ON:1
%    AcidosisFlag = 0        Acidosis OFF:0 / ON:1
%    HypoxiaIKatpFlag = 0    Hypoxia IKatp OFF:0 / ON:1
%    HypoxiaINaKFlag = 0     Hypoxia on INaK OFF:0 / ON:1
%    HypoxiaIpCaFlag = 0     Hypoxia on IpCa OFF:0 / ON:1
%    HypoxiaJupFlag = 0      Hypoxia on Jup OFF:0 / ON:1
%    HypoxiaCEFlag = 0       Hypoxia on CE OFF:0 / ON:1
%    IKatpType = 'Kazbanov'  Ledezma / Kazbanov
%
%   Internal factors that can be modified
%    gNa_factor = 1           gNa factor
%    g_f_factor = 1           g_f factor
%    g_CaL_factor = 1         g_CaL factor
%    g_to_factor = 1          g_to factor
%    g_Ks_factor = 1          g_Ks factor
%    g_Kr_factor = 1          g_Kr factor
%    g_K1_factor = 1          g_K1 factor
%    kNaCa_factor = 1         kNaCa factor
%    PNaK_factor = 1          PNaK factor
%    g_PCa_factor = 1         g_PCa factor
%    GNaLmax_factor = 1       GNaLmax factor
%    RyRtauadapt_factor = 1   RyRtauadapt factor
%    RyRtauact_factor = 1     RyRtauact factor
%    RyRtauinact_factor = 1   RyRtauinact factor
%    g_irel_max_factor = 1    g_irel_max factor
%    tau_m_factor = 1         tau_m factor
%    tau_h_j_factor = 1       tau_h_j factor
%    tau_d_factor = 1         tau_d factor
%    tau_f1_2_factor = 1      tau_f1_2 factor
%    RyRchalf_factor = 1      RyRchalf factor
%
%   classicVSoptic
%    classicVSoptic = 0  Temperature = 37°C: 0 / Temperature = 21°C: 1
%
%   Contractile type
%    contrflag = 1  isometric: 0 / active contraction: 1
%
%   Override file
%    override_file = ""  Allows override of anything before running (even other arguments!)
%
%   Initial value for differential equations
%    Vm_0 = -0.070        Initial Value for Vm
%    Ca_SR_0 = 0.32       Initial Value for Ca_SR_0
%    Cai_0 = 0.0002       Initial Value for Cai_0
%    g_0 = 0              Initial Value for g_0
%    d_0 = 0              Initial Value for d_0
%    f1_0 = 1             Initial Value for f1_0
%    f2_0 = 1             Initial Value for f2_0
%    fCa_0 = 1            Initial Value for fCa_0
%    Xr1_0 = 0            Initial Value for Xr1_0
%    Xr2_0 = 1            Initial Value for Xr2_0
%    Xs_0 = 0             Initial Value for Xs_0
%    h_0 = 0.75           Initial Value for h_0
%    j_0 = 0.75           Initial Value for j_0
%    m_0 = 0              Initial Value for m_0
%    Xf_0 = 0.1           Initial Value for Xf_0
%    q_0 = 1              Initial Value for q_0
%    r_0 = 0              Initial Value for r_0
%    Nai_0 = 9.2          Initial Value for Nai_0
%    m_L_0 = 0            Initial Value for m_L_0
%    h_L_0 = 0.75         Initial Value for h_L_0
%    RyRa_0 = 0.3         Initial Value for RyRa_0
%    RyRo_0 = 0.9         Initial Value for RyRo_0
%    RyRc_0 = 0.1         Initial Value for RyRc_0
%    Nxb_0 = 0.9997       Non-permissive fraction in overlap region
%    XBpreR_0 = 0.0001    Fraction of XBs in pre-rotated state (Stiffness generation)
%    XBpostR_0 = 0.0001   Fraction of XBS in post-rotated state (force generation)
%    x_XBpreR_0 = 0       Strain of XBs in the pre-rotated state
%    x_XBpostR_0 = 0.007  Strain of XBs in the post-rotated state
%    TropCaL_0 = 0.0144   Ca bound to the low affinity troponin site
%    TropCaH_0 = 0.1276   Ca bound to the high affinity troponin sit
%    IntegF_0 = 0         Initial Value for IntegF_0
%    SL_0 = 1.9           Sarcomere length
%    O2o_0 = 0.133        Initial Value for O2o_0
%
%  Example:
%    [result, time] = run_hiMCES(simTime = 10, result = 'Vm, Nai, Cai');
%    subplot(311)
%    plot(time, result(:,1)), xlabel("time"), ylabel("mV")
%    subplot(312)
%    plot(time, result(:,2)), xlabel("time"), ylabel("mM")
%    subplot(313)
%    plot(time, result(:,3)), xlabel("time"), ylabel("mM")
%
%
%
% Authors: Mohamadamin Forouzandehmehr (model) & Ossi Noita (This matlab-function)
%
% See also MasterCompute_hiMCES
	arguments
		args.result = 'Vm'      % Values stored into the result defined with comma separated list, in order
		args.override_file = "" % Allows override of any internal before running
		args.simTime = 800      % Simulation time
		args.stimFlag = 0				% Stimulation OFF:0 / ON:1
		args.tDrugApplication = 1400  % Time to start applying the drug
		args.dose_uM = 0.3      % Drug dosage, allowed values are 0.3 / 2 / 10
		args.INaLRedMed = 1     % INaLRedMed OFF:0 / ON:1
		args.IKsRedMed = 1      % IKsRedMed OFF:0 / ON:1
		args.INaCaRedMed = 1    % INaCaRedMed OFF:0 / ON:1
		args.IfRedMed = 1       % IfRedMed OFF:0 / ON:1
		args.HyperkalemiaFlag = 0  % Hyperkalemia OFF:0 / ON:1
		args.AcidosisFlag = 0      % Acidosis OFF:0 / ON:1
		args.HypoxiaIKatpFlag = 0  % Hypoxia IKatp OFF:0 / ON:1
		args.HypoxiaINaKFlag = 0   % Hypoxia on INaK OFF:0 / ON:1
		args.HypoxiaIpCaFlag = 0   % Hypoxia on IpCa OFF:0 / ON:1
		args.HypoxiaJupFlag = 0    % Hypoxia on Jup OFF:0 / ON:1
		args.HypoxiaCEFlag = 0     % Hypoxia on CE OFF:0 / ON:1
		args.IKatpType = 'Kazbanov'  % Ledezma or Kazbanov
		args.gNa_factor = 1          % latin hybercube idx 1
		args.g_f_factor = 1          % latin hybercube idx 2
		args.g_CaL_factor = 1        % latin hybercube idx 3
		args.g_to_factor = 1         % latin hybercube idx 4
		args.g_Ks_factor = 1         % latin hybercube idx 5
		args.g_Kr_factor = 1         % latin hybercube idx 6
		args.g_K1_factor = 1         % latin hybercube idx 7
		args.kNaCa_factor = 1        % latin hybercube idx 8
		args.PNaK_factor = 1         % latin hybercube idx 9
		args.g_PCa_factor = 1        % latin hybercube idx 10
		args.GNaLmax_factor = 1      % latin hybercube idx 11
		args.RyRtauadapt_factor = 1  % latin hybercube idx 12
		args.RyRtauact_factor = 1    % latin hybercube idx 13
		args.RyRtauinact_factor = 1  % latin hybercube idx 14
		args.g_irel_max_factor = 1   % latin hybercube idx 16
		args.tau_m_factor = 1        % latin hybercube idx 17
		args.tau_h_j_factor = 1      % latin hybercube idx 18
		args.tau_d_factor = 1        % latin hybercube idx 19
		args.tau_f1_2_factor = 1     % latin hybercube idx 20
		args.RyRchalf_factor = 1     % latin hybercube idx 21
		args.classicVSoptic = 0      % 0: Temperature = 37°C , 1: Temperature = 21°C
		args.contrflag = 1           % 0 = isometric  1 = active contraction
		args.ischemiaFlag = 0        % 0 = no ischemia, 1 = SEV1, 2 = SEV2
		args.Vm_0 = -0.070            % Initial Value for Vm
		args.Ca_SR_0 = 0.32           % Initial Value for Ca_SR_0
		args.Cai_0 = 0.0002           % Initial Value for Cai_0
		args.g_0 = 0                  % Initial Value for g_0
		args.d_0 = 0                  % Initial Value for d_0
		args.f1_0 = 1                 % Initial Value for f1_0
		args.f2_0 = 1                 % Initial Value for f2_0
		args.fCa_0 = 1                % Initial Value for fCa_0
		args.Xr1_0 = 0                % Initial Value for Xr1_0
		args.Xr2_0 = 1                % Initial Value for Xr2_0
		args.Xs_0 = 0                 % Initial Value for Xs_0
		args.h_0 = 0.75               % Initial Value for h_0
		args.j_0 = 0.75               % Initial Value for j_0
		args.m_0 = 0                  % Initial Value for m_0
		args.Xf_0 = 0.1               % Initial Value for Xf_0
		args.q_0 = 1                  % Initial Value for q_0
		args.r_0 = 0                  % Initial Value for r_0
		args.Nai_0 = 9.2              % Initial Value for Nai_0
		args.m_L_0 = 0                % Initial Value for m_L_0
		args.h_L_0 = 0.75             % Initial Value for h_L_0
		args.RyRa_0 = 0.3             % Initial Value for RyRa_0
		args.RyRo_0 = 0.9             % Initial Value for RyRo_0
		args.RyRc_0 = 0.1             % Initial Value for RyRc_0
		args.Nxb_0 = 0.9997           % Non-permissive fraction in overlap region
		args.XBpreR_0 = 0.0001        % Fraction of XBs in pre-rotated state (Stiffness generation)
		args.XBpostR_0 = 0.0001       % Fraction of XBS in post-rotated state (force generation)
		args.x_XBpreR_0 = 0           % Strain of XBs in the pre-rotated state
		args.x_XBpostR_0 = 0.007      % Strain of XBs in the post-rotated state
		args.TropCaL_0 = 0.0144       % Ca bound to the low affinity troponin site
		args.TropCaH_0 = 0.1276       % Ca bound to the high affinity troponin sit
		args.IntegF_0 = 0             % Initial Value for IntegF_0
		args.SL_0 = 1.9               % Sarcomere length
		args.O2o_0 = 0.133            % Initial Value for O2o_0
		args.MaxStep = 1e-3     % Longest allowed step for ode
		args.InitialStep = 2e-5 % Initial step for ode
	end

simTime =	args.simTime;
stimFlag =	args.stimFlag;
tDrugApplication =	args.tDrugApplication;
dose_uM =	args.dose_uM;
INaLRedMed =	args.INaLRedMed;
IKsRedMed =	args.IKsRedMed;
INaCaRedMed =	args.INaCaRedMed;
IfRedMed =	args.IfRedMed;
HyperkalemiaFlag =	args.HyperkalemiaFlag;
AcidosisFlag =	args.AcidosisFlag;
HypoxiaIKatpFlag =	args.HypoxiaIKatpFlag;
HypoxiaINaKFlag =	args.HypoxiaINaKFlag;
HypoxiaIpCaFlag =	args.HypoxiaIpCaFlag;
HypoxiaJupFlag =	args.HypoxiaJupFlag;
HypoxiaCEFlag =	args.HypoxiaCEFlag;
IKatpSelectionString =	args.IKatpType;
gNa_factor =	args.gNa_factor;
g_f_factor =	args.g_f_factor;
g_CaL_factor =	args.g_CaL_factor;
g_to_factor =	args.g_to_factor;
g_Ks_factor =	args.g_Ks_factor;
g_Kr_factor =	args.g_Kr_factor;
g_K1_factor =	args.g_K1_factor;
kNaCa_factor =	args.kNaCa_factor;
PNaK_factor =	args.PNaK_factor;
g_PCa_factor =	args.g_PCa_factor;
GNaLmax_factor =	args.GNaLmax_factor;
RyRtauadapt_factor =	args.RyRtauadapt_factor;
RyRtauact_factor =	args.RyRtauact_factor;
RyRtauinact_factor =	args.RyRtauinact_factor;
g_irel_max_factor =	args.g_irel_max_factor;
tau_m_factor =	args.tau_m_factor;
tau_h_j_factor =	args.tau_h_j_factor;
tau_d_factor =	args.tau_d_factor;
tau_f1_2_factor =	args.tau_f1_2_factor;
RyRchalf_factor =	args.RyRchalf_factor;
classicVSoptic =	args.classicVSoptic;
contrflag =	args.contrflag;
ischemiaFlag =	args.ischemiaFlag;
Vm_0 =	args.Vm_0;
Ca_SR_0 =	args.Ca_SR_0;
Cai_0 =	args.Cai_0;
g_0 =	args.g_0;
d_0 =	args.d_0;
f1_0 =	args.f1_0;
f2_0 =	args.f2_0;
fCa_0 =	args.fCa_0;
Xr1_0 =	args.Xr1_0;
Xr2_0 =	args.Xr2_0;
Xs_0 =	args.Xs_0;
h_0 =	args.h_0;
j_0 =	args.j_0;
m_0 =	args.m_0;
Xf_0 =	args.Xf_0;
q_0 =	args.q_0;
r_0 =	args.r_0;
Nai_0 =	args.Nai_0;
m_L_0 =	args.m_L_0;
h_L_0 =	args.h_L_0;
RyRa_0 =	args.RyRa_0;
RyRo_0 =	args.RyRo_0;
RyRc_0 =	args.RyRc_0;
Nxb_0 =	args.Nxb_0;
XBpreR_0 =	args.XBpreR_0;
XBpostR_0 =	args.XBpostR_0;
x_XBpreR_0 =	args.x_XBpreR_0;
x_XBpostR_0 =	args.x_XBpostR_0;
TropCaL_0 =	args.TropCaL_0;
TropCaH_0 =	args.TropCaH_0;
IntegF_0 =	args.IntegF_0;
SL_0 =	args.SL_0;
O2o_0 =	args.O2o_0;
MaxStep =	args.MaxStep;
InitialStep =	args.InitialStep;

if args.override_file ~= ""
	% Override any default value
	load(args.override_file)
end

options = odeset('MaxStep', MaxStep,'InitialStep', InitialStep);
%% Myofilament original initial state
Yn = [Vm_0, Ca_SR_0, Cai_0, g_0, d_0, f1_0, f2_0, fCa_0, Xr1_0, Xr2_0, Xs_0, h_0, j_0, m_0, Xf_0, q_0, r_0, Nai_0, m_L_0, h_L_0, RyRa_0, RyRo_0, RyRc_0, Nxb_0, XBpreR_0, XBpostR_0, x_XBpreR_0, x_XBpostR_0, TropCaL_0, TropCaH_0, IntegF_0, SL_0, O2o_0];

%% Levosimendan simulation
dose_uM = dose_uM;
levo_doses =    [0.3,      1,    2,        3,    5,   10,        15, 20];
fkon_to_doses = [1.1,      1.21, 1.32,     1.39, 1.57, 1.88,     2.2, 2.51];
flv_to_doses =  [1.002385,  NaN, 1.118538, NaN,  NaN,  1.843903, NaN, NaN];

search_nans = isnan(levo_doses + fkon_to_doses + flv_to_doses);
levo_idx = find(levo_doses == dose_uM);

if isempty(levo_idx) | search_nans(levo_idx)
  error(['Incorrect "dose_uM" does, use instead one of following: \n',...
		num2str(levo_doses(~search_nans))] )
end

fkon = fkon_to_doses(levo_idx);
flv =  flv_to_doses(levo_idx);

INaFIC50 = 85.714; % uM
INaFHill = 1;
INaFRedMed = 1./(((dose_uM)./(INaFIC50)).^INaFHill+1);

ICaLIC50 = 28.45; % uM
ICaLHill = 1;
ICaLRedMed = 1./(((dose_uM)./(ICaLIC50)).^ICaLHill+1);

IKrIC50 = 22; % uM
IKrHill = 0.68;
IKrRedMed = 1./(((dose_uM)./(IKrIC50)).^IKrHill+1);

%% ischemia setting
myFlags = [HyperkalemiaFlag AcidosisFlag  HypoxiaIKatpFlag  HypoxiaINaKFlag  HypoxiaIpCaFlag  HypoxiaJupFlag  HypoxiaCEFlag];
% z takes the extent of inhibition of INaK & IpCa and calculates
% corresponding normalized O2s
% fnp takes the level of oxygantion (0 to 1) and calculates the corresponding
% inhibition for INaK & IpCa

z = 1; fnp = 0;                   % Physoxia
if ischemiaFlag == 1
syms x
z = solve (1-rho(x) == 0.2, x);
z = double(z);
fnp = 1 - rho(z);

else if ischemiaFlag == 2
syms x
z = solve (1-rho(x) == 0.31, x);
z = double(z);
fnp = 1 - rho(z);

end
end

myLatinHypercubeVector = [gNa_factor, g_f_factor, g_CaL_factor, g_to_factor, g_Ks_factor, g_Kr_factor, g_K1_factor, kNaCa_factor, PNaK_factor, g_PCa_factor, GNaLmax_factor, RyRtauadapt_factor, RyRtauact_factor, RyRtauinact_factor, 1, g_irel_max_factor, tau_m_factor, tau_h_j_factor, tau_d_factor, tau_f1_2_factor, RyRchalf_factor];
%% Integration

if args.override_file ~= ""
	% Override any result from the default values (e.f. Ym)
	load(args.override_file)
end

[t,Yc] = ode15s(@hiMCES,[0 simTime],Yn, options, contrflag, myLatinHypercubeVector, stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, fkon, flv);
Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
d = Yc(:, 5);
f1 = Yc(:, 6);
f2 = Yc(:, 7);
fCa = Yc(:, 8);
Xr1 = Yc(:, 9);
Xr2 = Yc(:, 10);
Xs = Yc(:, 11);
h = Yc(:, 12);
j = Yc(:, 13);
m = Yc(:, 14);
Xf = Yc(:, 15);
q = Yc(:, 16);
r = Yc(:, 17);
Nai  = Yc(:,18);
m_L = Yc(:, 19);
h_L = Yc(:, 20);
RyRa = Yc(:, 21);
RyRo = Yc(:, 22);
RyRc = Yc(:, 23);
Nxb = Yc(:, 24);
XBpreR = Yc(:, 25);
XBpostR = Yc(:, 26);
x_XBpreR = Yc(:, 27);
x_XBpostR = Yc(:, 28);
TropCaL = Yc(:, 29);
TropCaH = Yc(:, 30);
IntegF = Yc(:, 31);
SL = Yc(:, 32);
O2o = Yc(:, 33);

full_data = zeros(size(Yc,1), 38);

for i= 1:size(Yc,1)
	[~, dati] = hiMCES(t(i), Yc(i,:), contrflag, myLatinHypercubeVector, stimFlag, classicVSoptic, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, fkon, flv);
	full_data(i, :) = dati;
end
INa       = full_data(:, 1);
If        = full_data(:, 2);
ICaL      = full_data(:, 3);
Ito       = full_data(:, 4);
IKs       = full_data(:, 5);
IKr       = full_data(:, 6);
IK1       = full_data(:, 7);
INaCa     = full_data(:, 8);
INaK      = full_data(:, 9);
IpCa      = full_data(:, 10);
IbNa      = full_data(:, 11);
IbCa      = full_data(:, 12);
Irel      = full_data(:, 13);
Iup       = full_data(:, 14);
Ileak     = full_data(:, 15);
Istim     = full_data(:, 16);
E_K       = full_data(:, 17);
E_Na      = full_data(:, 18);
INaL      = full_data(:, 19);
Fps       = full_data(:, 20);
AT        = full_data(:, 21);
ATPase    = full_data(:, 22);
Lsarc     = full_data(:, 23);
Velo      = full_data(:, 24);
JCB       = full_data(:, 25);
PSP       = full_data(:, 26);         % SERCA phosphorylation;
IKATP     = full_data(:, 27);
ronak     = full_data(:, 28);
O2s       = full_data(:, 29);         % Bath/source O2 concentration
O2e       = full_data(:, 30);         % Extracellular O2 concentration
vcy       = full_data(:, 31);         % SERCA cycle rate
do2dt     = full_data(:, 32);
fhib      = full_data(:, 33);        % Returns finhib
Koo       = full_data(:, 34);
ro1       = full_data(:, 35);
ro2       = full_data(:, 36);
fkatp     = full_data(:, 37);
kon       = full_data(:, 38);

mat_currents = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IpCa, IbNa, IbCa, Irel, Iup, Ileak, Istim, INaL];
I_tot = sum(mat_currents, 2);

if isstring(args.result)
	args.result = convertStringsToChars(args.result);
end
result = eval(['[', args.result, ']']);
% End
