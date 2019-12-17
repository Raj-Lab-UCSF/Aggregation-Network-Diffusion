function franchi_raj_cleaned

% Code for the Aggregation-Network Diffusion (AND) model of pathology
% ramification in human brain connectome. Please cite the paper:
% Combined Model of Aggregation And Network Diffusion Recapitulates Alzheimer’s Regional Tau-PET 
% Ashish Raj*, Veronica Tora, Xiao Gao, Hanna Cho, Jae Yong Choi, Young Hoon Ryu, Chul Hyoung Lyoo, Bruno Franchi
% *Department of Radiology and Biomedical Imaging, University of California at San Francisco

% This version contains code that can do both tau and Abeta
% For now we are enabling only dynamics of tau, no Ab (will add later for paper 2)
% This code will reproduce the result sof the fully realised Franchi-Raj AND model
% We are posting group regional atrophy and tau SUVr tables for use in this
% code, however we are unable to post individual subject data

% Please add all subfolders within this folder to matlab path prior to running this code

suppress_glassbrain = 0; % use this to suppress the slow glassbrain rendering process (for testing/debugging)
orig_NDM_flag = 0;  % Set to 0 for AND model, 1 for original NDM (cannot do both...)
figuresavename = 'test1';
sig = 2;
sig_abeta = 2;
rescale_method = 'logistic';
ballradius = 5;  % just for glassbrain plotting
node_scrambling = 0;  % permutation testing of ODE - set to 1 if needed
edge_scrambling = 0;  % permutation testing on random connectomes - set to 1 if needed
corrtype = 'Pearson'; % 'Pearson' 'Spearman' or 'Overlap' or 'Dotproduct' or 'boxcox'; was DotProduct
pthr = 0.5; % only used for corrtype overlap; higher is more AD specific

fontsz = 16; % contrl sfont size on figs and plots
% Atlas, connectome, etc - incomplete revision, do later!
which_atlas = 'DK86'; % 'FXCN' or 'DK86'

thisdir = pwd;
permHCP = [19:52, 53:86, 1:9, 10:18];  % permutes 86 x 86 HCP conn matrix so that subcorts are at the end, as in previous work
    Cdk_hcp = csvread('mean80_fibercount.csv', 1, 0);         
    Cdk_hcp = Cdk_hcp(permHCP, permHCP);
    C = Cdk_hcp;
nroi = size(C,1);

% Define some useful indices
ECLind = 5; ECRind = 34+ECLind;  
AMLind = 75; AMRind = 84; 
CDLind = 71; CDRind = 80; 
PULind = 72; PURind = 81;
latOCCLind = 10; latOCCRind = 44;
ThLind = 70; ThRind = 79;
% Striatal_inds = [70 71 72 73 76 79 80 81 82 85];
Striatal_inds = [71 72 73 76 80 81 82 85];  % included thalamus
cerebrum_rois = 1:86;
noncerebrum_rois = 87:90;
ii = setdiff(cerebrum_rois, Striatal_inds);  % all regions not in striatum or cerebellum

s = load ('FSatlas_lobes.mat');
lobes = s.lobes;

% condition the matrix
thr = 5*mean(C(C>0));
C = min(C, thr);

% spectrum stuff
% L = diag(rowdegree) - C;
% L = eye(nroi) - diag(1./(rowdegree+eps)) * C;
% L = eye(nroi) - diag(1./(sqrt(rowdegree)+eps)) * C* diag(1./(sqrt(coldegree)+eps)) ;
% L = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * C; % * diag(1./(sqrt(coldegree)+eps)) ;
% L = eye(nroi) - diag(1./(sqrt(sqrt(rowdegree.*coldegree)+eps))) * C * diag(1./(sqrt(sqrt(rowdegree.*coldegree)+eps)));
% L = eye(nroi) - diag(1./(sqrt((rowdegree+coldegree)/2)+eps)) * C; % * diag(1./(sqrt(coldegree)+eps)) ;


% Do both retro and antero
Cant = C.'; %antero
% Cant = (C+C.')/2; %bi, remove latr
rowdegree = (sum(Cant, 2)).';
coldegree = sum(Cant, 1);
Lant = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Cant; % * diag(1./(sqrt(coldegree)+eps)) ;
% Lant = eye(nroi) - diag(1./(sqrt(rowdegree)+eps)) * Cant* diag(1./(sqrt(coldegree)+eps)) ;
[V, D] = eig(Lant); 
[dd, id] = sort(diag(abs(D)), 'Ascend');
evalues_ant = dd;
evecs_ant = V(:,id);

Cret = C; %retro
rowdegree = (sum(Cret, 2)).';
coldegree = sum(Cret, 1);
Lret = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Cret; % * diag(1./(sqrt(coldegree)+eps)) ;
% Lret = eye(nroi) - diag(1./(sqrt(rowdegree)+eps)) * Cret * diag(1./(sqrt(coldegree)+eps)) ;
[V, D] = eig(Lret); 
[dd, id] = sort(diag(abs(D)), 'Ascend');
evalues_ret = dd;
evecs_ret = V(:,id);


%% Load some real data

% APP map, artificially created from literature 
% (occ and parital are 0.5 of others, striatals are 1.1x)
nc_app = ones(nroi, 1);
nc_app(lobes==3) = 0.5;
nc_app(lobes==1) = 1.1;
nc_app([CDLind, CDRind, PULind, PURind]) = 1.1;
% MAPT map, artificially created from literature (occ and striatal are 0.5 of others)
nc_mapt = ones(nroi, 1);
nc_mapt(lobes==3) = 0.3;
nc_mapt([CDLind, CDRind, PULind, PURind]) = 0.5;

% load ADNI group stats
% cd sneha_adni_group;
cd Fon_ADNI2_Desikan_Group_Stats

S = load('Age_matched_stats.mat');
nc_met = S.CON1_FDG_baseline(:,3);
nc_met(isnan(nc_met)) = 0;
emci_met = S.EMCI_baseline_FDG_stats(:,3);
emci_met(isnan(emci_met)) = 0;
lmci_met = S.LMCI_baseline_FDG_stats(:,3);
lmci_met(isnan(lmci_met)) = 0;
ad_met = S.AD_baseline_FDG_stats(:,3);
ad_met(isnan(ad_met)) = 0;

nc_abeta = S.CON1_AV45_baseline_stats(:,3);
nc_abeta(isnan(nc_abeta)) = 0;
emci_abeta =  -S.EMCI_baseline_AV45_stats(:,3);
emci_abeta(isnan(emci_abeta)) = 0;
lmci_abeta = -S.LMCI_baseline_AV45_stats(:,3);
lmci_abeta(isnan(lmci_abeta)) = 0;
ad_abeta = -S.AD_baseline_AV45_stats(:,3);
ad_abeta(isnan(ad_abeta)) = 0;

nc_atrophy = -S.CON1_vol_1_stats(:,3);
nc_atrophy(isnan(nc_atrophy)) = 0;

emci_atrophy = S.EMCI_vol_1_stats(:,3);
emci_atrophy(isnan(emci_atrophy)) = 0;
lmci_atrophy = S.LMCI_vol_4_stats(:,3);
lmci_atrophy(isnan(lmci_atrophy)) = 0;

ad_atrophy = S.AD_vol_1_stats(:,3);
ad_atrophy(isnan(ad_atrophy)) = 0;

% Collect amyloid +/- subjects in separate groups for specific hypothesis
% testing regarding amyloid status
S = load('Pos_Neg_FinalStats.mat');
mciabpos_met = -S.stats_POS_CON_FDG(:,3);
mciabpos_atrophy = -S.stats_POS_CON_Vol(:,3);
mciabpos_abeta = S.stats_POS_CON_AV45(:,3);
mciabneg_met = -S.stats_NEG_CON_FDG(:,3);
mciabneg_atrophy = -S.stats_NEG_CON_Vol(:,3);
mciabneg_abeta = S.stats_NEG_CON_AV45(:,3);
mciabpos_neg_met = -S.stats_POS_NEG_FDG(:,3);
mciabpos_neg_atrophy = -S.stats_POS_NEG_Vol(:,3);
mciabpos_neg_abeta = S.stats_POS_NEG_AV45(:,3);

S = load('Tau_Structure.mat');
taudata = S.Tau_Structure;
patient_tau = [];
for i = 1:length(taudata)
    q = taudata(i);
    disp(q.DX);
%     if strcmp(q.DX, 'EMCI') || strcmp(q.DX, 'LMCI') || strcmp(q.DX, 'MCI') || strcmp(q.DX, 'AD') 
    if ~strcmp(q.DX, 'Normal')
        patient_tau = [patient_tau, q.Tau(:)];
    end  
end
patient_tau = mean(patient_tau - 1,2);

cd(thisdir);

% Import Korea tau data from Yonsei U
% Only group stats regional data are being shared here; this should be adequate for the publishec paper
% for individual subjects please inquire directly with Yonsei University researchers

qq = xlsread('Group_T807_Tau.xlsx');

q = qq(4,2:end);
korea_namci_tau = q(:) - 1;
q = qq(3,2:end);
korea_amci_tau = q(:) - 1;
q = qq(2,2:end);
korea_ad_tau = q(:) - 1;

qq = xlsread('Group_NRC_Amyloid.xlsx');
q = qq(4,2:end);
korea_namci_abeta = q(:) - 1;
q = qq(3,2:end);
korea_amci_abeta = q(:) - 1;
q = qq(2,2:end);
korea_ad_abeta = q(:) - 1;

% rescale brain data to convert real axis of atrophy/tau etc to [0,1]
% We have implemenbetd a logistic transform function for this, controlled
% by the width parameter sig
temp = rescale_braindata([nc_met emci_met lmci_met ad_met], sig, rescale_method, 'allcols');
nc_met = temp(:,1);
emci_met = temp(:,2);
lmci_met = temp(:,3);
ad_met = temp(:,4);
temp = rescale_braindata([nc_abeta emci_abeta lmci_abeta ad_abeta], sig_abeta, rescale_method, 'allcols');
nc_abeta = temp(:,1);
emci_abeta = temp(:,2);
lmci_abeta = temp(:,3);
ad_abeta = temp(:,4);
temp = rescale_braindata([emci_atrophy lmci_atrophy ad_atrophy], sig, rescale_method, 'allcols');
emci_atrophy = temp(:,1);
lmci_atrophy = temp(:,2);
ad_atrophy = temp(:,3);
temp = rescale_braindata([mciabpos_met mciabneg_met mciabpos_neg_met], sig, rescale_method, 'allcols');
mciabpos_met = temp(:,1);
mciabneg_met = temp(:,2);
mciabpos_neg_met = temp(:,3);
temp = rescale_braindata([mciabpos_abeta mciabneg_abeta mciabpos_neg_abeta], sig, rescale_method, 'allcols');
mciabpos_abeta = temp(:,1);
mciabneg_abeta = temp(:,2);
mciabpos_neg_abeta = temp(:,3);
temp = rescale_braindata([mciabpos_atrophy mciabneg_atrophy mciabpos_neg_atrophy], sig, rescale_method, 'allcols');
mciabpos_atrophy = temp(:,1);
mciabneg_atrophy = temp(:,2);
mciabpos_neg_atrophy = temp(:,3);

temp = rescale_braindata([korea_namci_abeta korea_amci_abeta korea_ad_abeta], sig_abeta, rescale_method, 'allcols');
korea_namci_abeta = temp(:,1);
korea_amci_abeta = temp(:,2);
korea_ad_abeta = temp(:,3);
% figure; plot(1:86, korea_namci_abeta, 1:86, korea_amci_abeta, 1:86, korea_ad_abeta),

% figure; plot(1:86, korea_namci_tau, 1:86, korea_amci_tau, 1:86, korea_ad_tau),
temp = rescale_braindata([korea_namci_tau korea_amci_tau korea_ad_tau], sig, rescale_method, 'allcols');
korea_namci_tau = temp(:,1);
korea_amci_tau = temp(:,2);
korea_ad_tau = temp(:,3);
% figure; plot(1:86, korea_namci_tau, 1:86, korea_amci_tau, 1:86, korea_ad_tau),

patient_tau = rescale_braindata(patient_tau, sig, rescale_method, 'allcols');


%% 7) Raj_Franchi AND model (Aggregation and Network Diffusion)

    % Implements the Raj-Franchi model of agglomeration and network diffusion
    M = 5; % number of aligomers, at M tangles/placques form
    sM = 2; %3; % a small number, width of the aggregation reaction curve
%     A = ones(M,M);  % coeffs a_ij, change later
%     A = M./((1:M).'*(1:M));  % coeffs a_ij= 1/ij, based on statistical mechanics
        A = ((1:M).'*(1:M))/sM/sM .* exp(-(1:M).'*(1:M)/sM/sM);  % coeffs a_ij= gamma(ij) function, based on empirical observations
    dconst = 1; %4/M;   % parameter of exponential decay mode governing the diffusion constants of oligomers of varying length
    
    % diffvec is the relationship between diffusivity and oligomer size.
    % This could be any monotonically decreasing function:
    diffvec = exp(-dconst*(1:M).'); % exp decay
%     diffvec =  2*((1:M).^(-2)).'; % inverse power law
%     diffvec = (1:M).'/sM .* exp(-(1:M).'/sM);  % gamm function

    tmax = 400;
    trange = [0,tmax];
    tsamples = [5 10 15 20]; %[5 10 15 25];
%     xdrive = 400*abs(evecs_ret(:,2)).*nc_mapt;  %10
%     xdrive ( xdrive < mean(xdrive(:)) + std(xdrive(:)) ) = 0;  
%     xdrive ( xdrive < mean(xdrive(:)) ) = 0;  

% Constantrs below are optimized via exhaustive search (one time run).
% Below are their optimum values
    cAb = 10; % global constant weighing Ab diffusion against production and agglomeration terms
    ctau = 4; %20; % global constant weighing tau diffusion against production and agglomeration terms       
    beta = 1; %20; %4
    alpha = 0.05; %1*beta; %0.5
    % gama = 8*beta; %0.8 8
    gama = 0; % enforce no interaction
    lambda = 15;
    
% Define the seed location below. 
    xdrive = zeros(nroi,1);
    xdrive([ECLind, ECRind]) = 800;  % EC seeding
%     xdrive([nroi-1, nroi]) = 500;  % LC seeding
%     xdrive([PULind, PURind]) = 200;  % EC seeding
    % Note, 800 is a magic number, needed simply for the brain plotting
    % routines. There is no inherent reason for it to exist
    
    % define the driving term that controls the production of Abeta as a function of baseline
    % metabolic activity of each region, as well as the local pool of available APP
    fdrive = 0.5*nc_met.*nc_app;
    x0 = zeros([size(xdrive), M]);
    f0 = zeros([size(fdrive), M]);
    xf0 = [ x0(:); f0(:)];  
    
if orig_NDM_flag == 0
    % Now solve the ODE using non-negative constraints in the time range[0, trange]
    % Note: We have set up a joint model of tau and Abeta, with a future
    % project in mind where the two interact. For the PLoS CB paper the
    % interaction term is disabled and teh Abeta variables are ignored
    odeopts = odeset('NonNegative', 1:length(xf0), 'RelTol',1e-2,'AbsTol',1e-4);
    [tsol, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
    nt = size(xfsol_all,1);
    xsol_all = xfsol_all(:, 1:nroi*M); 
    xsol_all = reshape(xsol_all, [nt, nroi,M]);
    xsol_all = permute(xsol_all, [2,1,3]);  
    xsol = xsol_all(:,:,M);
    fsol_all = xfsol_all(:, nroi*M+1:2*nroi*M); 
    fsol_all = reshape(fsol_all, [nt, nroi,M]);
    fsol_all = permute(fsol_all, [2,1,3]);
    fsol = fsol_all(:,:,M);

    % comment out as not using interaction terms
%     gama = 0;
%     [tsolnoint, xfsolnoint_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%     nt = size(xfsolnoint_all,1);
%     xsolnoint_all = xfsolnoint_all(:, 1:nroi*M); 
%     xsolnoint_all = reshape(xsolnoint_all, [nt, nroi,M]);
%     xsolnoint_all = permute(xsolnoint_all, [2,1,3]);    
%     xsolnoint = xsolnoint_all(:,:,M);
%     fsolnoint_all = xfsolnoint_all(:, nroi*M+1:2*nroi*M); 
%     fsolnoint_all = reshape(fsolnoint_all, [nt, nroi,M]);
%     fsolnoint_all = permute(fsolnoint_all, [2,1,3]);
%     fsolnoint = fsolnoint_all(:,:,M);

elseif orig_NDM_flag == 1
% DEBUG: run original NDM, with no aggregation
    tsol = trange(1):trange(end);
    nt = length(tsol);
    bet = beta*0.2;
    xsol_all = zeros(nroi,nt, M);
    xsol = zeros(nroi,nt);
    for qi = 1:nt
        q = 0.05*sqrt(tsol(qi))*expm(-Lret*bet*tsol(qi))*xdrive;  % magic multiple 10 needed to ensure glassbrain balls the right size... ugly but works
        % Note: sqrt above is purely for visualization; to keep early
        % snapshots from giving too-high ball sizes. This is expected since
        % NDM is mass preserving
        xsol(:,qi) = q;
        xsol_all(:,qi,:) = q*ones(1,M);
    end
    fsol_all = xsol_all;
    fsol = xsol;
end
save xsol.mat xsol;

     % Global pathology burden
      for k = 1:M
        globaltau(k,:) = sum(xsol_all(:,:,k),1);
        globalab(k,:) = sum(fsol_all(:,:,k),1);
      end
     globaltau = globaltau.'; 
     globalab = globalab.';
     CSFtau = globaltau*diffvec;  % added 2/19
     CSFab = globalab*diffvec;    % added 2/19
     
    % % Plot Rt curves w.r.t. ADNI volumetric data
    figure; 
    subplot(121); plot(tsol, globalab);   title('Global burden of modeled A\beta');
    Rnc = mycorr(fsol(cerebrum_rois,:), nc_abeta(cerebrum_rois));
    Remci = mycorr(fsol(cerebrum_rois,:), emci_abeta(cerebrum_rois));
    Rlmci = mycorr(fsol(cerebrum_rois,:), lmci_abeta(cerebrum_rois));
    Rad = mycorr(fsol(cerebrum_rois,:), ad_abeta(cerebrum_rois));
    subplot(122); plot(tsol, Rnc, 'b-', tsol, Remci, 'r:', tsol, Rlmci, 'r-', tsol, Rad, 'k-'); title('R, empirical A\beta against ODE solution');

    h = plot(tsol, [globaltau, CSFtau]);  legend('m=1 (monomer)', 'm=2', 'm=3', 'm=4', 'm=5 (tangle)', 'CSF tau');
    [h.LineWidth] = deal(2); h(6).LineWidth = 3;
    hold on; 
    title('Global burden of modeled tau', 'FontSize', fontsz); xlabel('model time (a.u.)', 'FontSize', fontsz); ylabel('\tau_m(t)', 'FontSize', fontsz);
    Ramci = [];
    for m = 1 :M 
       Ramci(m,:) = mycorr(xsol_all(ii,:,m), korea_amci_tau(ii)); 
    end
    subplot(132); 
    h = plot(tsol, Ramci); [h.LineWidth] = deal(3); 
    title('R, AV1451-PET of aMCI against ODE-predicted oligomers', 'FontSize', fontsz); xlabel('model time (a.u.)', 'FontSize', fontsz); ylabel('R', 'FontSize', fontsz);
    Rlmci = [];
    for m = 1 :M 
       Rlmci(m,:) = mycorr(xsol_all(cerebrum_rois,:,m), lmci_atrophy(cerebrum_rois)); 
    end
    subplot(133); 
    h = plot(tsol, Rlmci); 
    [h.LineWidth] = deal(3); 
    title('R, atrophy against ODE solution', 'FontSize', fontsz); xlabel('model time (a.u.)', 'FontSize', fontsz); ylabel('R', 'FontSize', fontsz);

    
%% Added 2/19 - fit CSF tau to data
    
    % Table below was created using group average data of CSF biomarkers from ADNI
    % Manually tweaked scale constants so that the scale of AND will match
    % teh scale of teh CSF-tau raw values. 
    % Similarly, the time of each event (e.g. aMCI group) is unknown, and
    % was tehrefore manually tweaked within a reasonable range. The purpose
    % here is not rigorous model fitting, but simply to show that at a
    % "good enough" scale and shift value the AND-derived biomarker begins
    % to replicate teh CSF-tau emirical data
    
    ADNI_CSFptau = [21.149, 24.319, 30.514, 36.162];  % [sSMC, EMCI, LMCI, AD]
    ADNI_CSFptau_std = [9.993 13.394 14.784 15.115];  
    ADNI_CSFptau_slope = [0.725 0.320 0.323 -0.733];
     % hand-drawn model time for each Dx - this was tweaked manually using
     % trial and error, until the best fit with AND model was achieved
%     modeltime_perDx = tmax*[0.17, 0.22, 0.30, 0.7]; 
    modeltime_perDx = tmax*[0.12, 0.16, 0.22, 0.50];  % hand-drawn model time for each Dx
    
    % Now add converters v non-cnverters
    ADNI_CSFptau_conv = [25.365 28.503 35.423 34.825];  % [pHC, pSMC, pEMCI, pLMCI]
    ADNI_CSFptau_std_conv = [9.401 14.358 19.207 14.760];  
    ADNI_CSFptau_slope_conv = [0.725 1.817 0.296 0.320];
    modeltime_perDx_conv(1) = (modeltime_perDx(1) + modeltime_perDx(2))/2;  % place time mid-way between Dx times
    modeltime_perDx_conv(2) = (modeltime_perDx(1) + modeltime_perDx(3))/2;  % place time mid-way between Dx times
    modeltime_perDx_conv(3) = (modeltime_perDx(2) + modeltime_perDx(4))/2;  % place time mid-way between Dx times
    modeltime_perDx_conv(4) = (modeltime_perDx(3) + modeltime_perDx(4))/2;  % place time mid-way between Dx times
 
    % Fit adni = const*CSFtau to get the best global scale factor
    for kq = 1:length(ADNI_CSFptau)
       [~,ind] = min(abs(tsol-modeltime_perDx(kq)));
       qtemp(kq) = CSFtau(ind);       
    end
    for kq = 1:length(ADNI_CSFptau_conv)
       [~,ind] = min(abs(tsol-modeltime_perDx_conv(kq)));
       qtemp(kq+length(ADNI_CSFptau)) = CSFtau(ind);       
    end
    newCSFtau = norm([ADNI_CSFptau ADNI_CSFptau_conv])/norm(qtemp)*CSFtau;

%     % Fit a polynomial p of degree 1 to the (x,y) data: discontinued
%     p = polyfit(qtemp,ADNI_CSFptau,1);
%     % Evaluate the fitted polynomial p and plot:
%     newCSFtau = polyval(p,CSFtau);
    
figure; hold on
plot(tsol, newCSFtau, 'b-','LineWidth', 3); 
errorbar(modeltime_perDx, ADNI_CSFptau, ADNI_CSFptau_std, 's','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red'); 
errorbar(modeltime_perDx_conv, ADNI_CSFptau_conv, ADNI_CSFptau_std_conv, 's','MarkerSize',10, 'MarkerEdgeColor','red'); 

% scatter(modeltime_perDx, ADNI_CSFptau, 100, 'r', 's','filled'); 
% scatter(modeltime_perDx_conv, ADNI_CSFptau_conv, 100, 'r', 's'); 
tseg = linspace(-2,2, 10);
for kq = 1:length(ADNI_CSFptau)
    t0 = modeltime_perDx(kq);
    y0 = ADNI_CSFptau(kq);
    sl = ADNI_CSFptau_slope(kq);
    yseg = sl*(tseg)+y0;
    plot(t0+tseg, yseg, 'k-', 'LineWidth', 2);    
end
for kq = 1:length(ADNI_CSFptau_conv)
    t0 = modeltime_perDx_conv(kq);
    y0 = ADNI_CSFptau_conv(kq);
    sl = ADNI_CSFptau_slope_conv(kq);
    yseg = sl*(tseg)+y0;
    plot(t0+tseg, yseg, 'k-', 'LineWidth', 2);    
end

hold off;

% Now add individual ADNI subjects CSF data
CSF_filename = 'CSF_Protein_regression-axis.xlsx';
    CSFdata = xlsread(CSF_filename, 'Data_Summation');

nsubj = size(CSFdata,1);
pid = CSFdata(:,1);
ADNI_CSFptau = CSFdata(:,7);
ADNI_CSFptau_z = CSFdata(:,8);
ADNI_CSF_Dx = CSFdata(:,2);
% remove longitudinal
[cc, ind] = unique(pid, 'last');
pid = pid(ind);
nsubj = length(pid);
ADNI_CSFptau = ADNI_CSFptau(ind);
ADNI_CSFptau_z = ADNI_CSFptau_z(ind);
ADNI_CSF_Dx = ADNI_CSF_Dx(ind);

z_shift = 7; % shift empirical z-score to match the model time at origin - manually tweaked until best performance achieved using trial and error. This may be replaced by a more principled optimization procedure
z_scale = 0.6; % scale empirical z-score to match the model time range - also manually optimized
ADNI_CSFptau_z= z_scale*ADNI_CSFptau_z + z_shift ;

% Fit adni = const*CSFtau    
qtemp = zeros(nsubj,1);
for kq = 1:nsubj
   thisz = ADNI_CSFptau_z(kq);
   [~,ind] = min(abs(tsol-thisz));
   qtemp(kq) = CSFtau(ind);       
end
newCSFtau = norm(ADNI_CSFptau)/norm(qtemp)*CSFtau;
qtemp = norm(ADNI_CSFptau)/norm(qtemp)*qtemp;
[R,p] = corr(ADNI_CSFptau, qtemp);

figure; hold on; 
sz = 40;
scatter(ADNI_CSFptau_z, ADNI_CSFptau, sz, ADNI_CSF_Dx, 'o', 'filled'); 
colormap(jet(5)); colorbar;
plot(tsol, newCSFtau, 'b-','LineWidth', 4);  
hold off;


%% % 6-7a) Plot Rt curves w.r.t. ADNI and Korea regional data

% For now, commenting out Ab results
% % First the ODE dynamics of fsol and empirical Abeta
%     
% for i = 1:length(tsamples)
%     figstr = [figuresavename '-ODEf' num2str(i)];
%     [dummy, ind] = min(abs(tsol-tsamples(i)));
%     yy = 0.5*ballradius*fsol(:,ind);  %5
% %     yy = min(yy, mean(yy) + 1*std(yy));
%     %yy = max(0, yy - 1*std(yy));
%     [h, tiffnamef{i}] = glassbrain_display(yy(cerebrum_rois), figstr); 
% end
% [Rall_emci, R_emci] = mycorr(fsol(cerebrum_rois,:), emci_abeta(cerebrum_rois));
% % [~, ind_emci] = max(Rall_emci);
% % [R_emci] = mycorr(fsol(cerebrum_rois, ind_emci), abs(emci_abeta(cerebrum_rois)));
% [Rall_lmci, R_lmci] = mycorr(fsol(cerebrum_rois,:), lmci_abeta(cerebrum_rois));
% % [~, ind_lmci] = max(Rall_lmci);
% % [R_lmci] = mycorr(fsol(cerebrum_rois, ind_lmci), abs(lmci_abeta(cerebrum_rois)));
% [Rall_ad, R_ad] = mycorr(fsol(cerebrum_rois,:), ad_abeta(cerebrum_rois));
% % [~, ind_ad] = max(Rall_ad);
% % [R_ad] = mycorr(fsol(cerebrum_rois, ind_ad), abs(ad_abeta(cerebrum_rois)));
% 
% [Rall_korea_namci, R_korea_namci] = mycorr(fsol(cerebrum_rois,:), korea_namci_abeta(cerebrum_rois));
% % [~, ind_korea_namci] = max(Rall_korea_namci);
% % [R_korea_namci] = mycorr(fsol(cerebrum_rois, ind_korea_namci), abs(korea_namci_abeta(cerebrum_rois)));
% [Rall_korea_amci, R_korea_amci] = mycorr(fsol(cerebrum_rois,:), korea_amci_abeta(cerebrum_rois));
% % [~, ind_korea_amci] = max(Rall_korea_amci);
% % [R_korea_amci] = mycorr(fsol(cerebrum_rois, ind_korea_amci), abs(korea_amci_abeta(cerebrum_rois)));
% [Rall_korea_ad, R_korea_ad] = mycorr(fsol(cerebrum_rois,:), korea_ad_abeta(cerebrum_rois));
% % [~, ind_korea_ad] = max(Rall_korea_ad);
% % [R_korea_ad] = mycorr(fsol(cerebrum_rois, ind_korea_ad), abs(korea_ad_abeta(cerebrum_rois)));
% 
%         figstr = [figuresavename '-emciab'];
%     yy =5*ballradius*emci_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname2] = glassbrain_display(yy(cerebrum_rois), figstr);
%         figstr = [figuresavename '-korea_namciab'];
%     yy =ballradius*korea_namci_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname3] = glassbrain_display(yy(cerebrum_rois), figstr);
%         figstr = [figuresavename '-lmciab'];
%     yy =5*ballradius*lmci_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname4] = glassbrain_display(yy(cerebrum_rois), figstr);
%         figstr = [figuresavename '-korea_amciab'];
%     yy =ballradius*korea_amci_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname5] = glassbrain_display(yy(cerebrum_rois), figstr);
%         figstr = [figuresavename '-adab'];
%     yy =5*ballradius*ad_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname6] = glassbrain_display(yy(cerebrum_rois), figstr);
%         figstr = [figuresavename '-korea_adab'];
%     yy =ballradius*korea_ad_abeta;
%     yy(Striatal_inds) = min(yy);
%     [h, tiffname7] = glassbrain_display(yy(cerebrum_rois), figstr);
% 
% figure; 
% cnt = 0;
% for i = 1:length(tsamples)
%     subplot(length(tsamples),3,cnt+1); subimage(imread(tiffnamef{i}));  title(['A\beta spread, t = ' num2str(tsamples(i))]); axis off
%     cnt = cnt+3;
% end
% % skiprows = length(tsamples) - 4;
% % cnt = 3*skiprows;
% % subplot(length(tsamples),3,cnt+2); scatter(fsol(cerebrum_rois, ind_emci), abs(emci_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_emci)]); xlabel('Model at peak similarity'); ylabel('ADNI EMCI AV45 uptake');
% % subplot(length(tsamples),3,cnt+3); scatter(fsol(cerebrum_rois, ind_korea_namci), abs(korea_namci_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_korea_namci)]); xlabel('Model at peak similarity'); ylabel('Korea naMCI AV45 uptake');
% % cnt = cnt + 3;
% % subplot(length(tsamples),3,cnt+2); scatter(fsol(cerebrum_rois, ind_lmci), abs(lmci_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_lmci)]); xlabel('Model at peak similarity'); ylabel('ADNI LMCI AV45 uptake');
% % subplot(length(tsamples),3,cnt+3); scatter(fsol(cerebrum_rois, ind_korea_amci), abs(korea_amci_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_korea_amci)]); xlabel('Model at peak similarity'); ylabel('Korea aMCI AV45 uptake');
% % cnt = cnt + 3;
% % subplot(length(tsamples),3,cnt+2); scatter(fsol(cerebrum_rois, ind_ad), abs(ad_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_ad)]); xlabel('Model at peak similarity'); ylabel('ADNI AD AV45 uptake');
% % subplot(length(tsamples),3,cnt+3); scatter(fsol(cerebrum_rois, ind_korea_ad), abs(korea_ad_abeta(cerebrum_rois)), 50, 'b', 'fill'); lsline; 
% % title(['R = ' num2str(R_korea_ad)]); xlabel('Model at peak similarity'); ylabel('Korea AD AV45 uptake');
% 
% cnt = 0;
% subplot(length(tsamples),3,cnt+2); plot(tsol, Rall_emci, 'g-', tsol, Rall_lmci, 'r-', tsol, Rall_ad, 'b-');
% title('Similarity index between model and ADNI patients AV45-PET (green: EMCI, red: LMCI, blue: AD)'); xlabel('Model time (a.u.)'); ylabel('Similarity index');
% subplot(length(tsamples),3,cnt+3); plot(tsol, Rall_korea_namci, 'g-', tsol, Rall_korea_amci, 'r-', tsol, Rall_korea_ad, 'b-');
% title('Similarity index between model and Korea patients AV45-PET (green: naMCI, red: aMCI, blue: AD)'); xlabel('Model time (a.u.)'); ylabel('Similarity index');
% cnt = 3;
% subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname2));  title(['ADNI EMCI AV45 uptake, Rmax = ' num2str(R_emci)]); axis off
% subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname3));  title(['Korea naMCI AV45 uptake, Rmax = ' num2str(R_korea_namci)]); axis off
% cnt = cnt+3;
% subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname4));  title(['ADNI LMCI AV45 uptake, Rmax = ' num2str(R_lmci)]); axis off
% subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname5));  title(['Korea aMCI AV45 uptake, Rmax = ' num2str(R_korea_amci)]); axis off
% cnt = cnt + 3;
% subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname6));  title(['ADNI AD AV45 uptake, Rmax = ' num2str(R_ad)]); axis off
% subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname7));  title(['Korea AD AV45 uptake, Rmax = ' num2str(R_korea_ad)]); axis off
 

% Second, ODE dynamics of xsol and empirical atrophy, tau
for i = 1:length(tsamples)
    figstr = [figuresavename '-ODEx' num2str(i)];
    [dummy, ind] = min(abs(tsol-tsamples(i)));
    yy = 0.2*ballradius*xsol(:,ind); % 0.1
%     yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamex{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
end
[Rall_emci, R_emci] = mycorr(xsol(cerebrum_rois,:), emci_atrophy(cerebrum_rois));
[~, ind_emci] = max(Rall_emci);
% [R_emci] = mycorr(xsol(cerebrum_rois, ind_emci), abs(emci_atrophy(cerebrum_rois)));
[Rall_lmci, R_lmci] = mycorr(xsol(cerebrum_rois,:), lmci_atrophy(cerebrum_rois));
[~, ind_lmci] = max(Rall_lmci);
% [R_lmci] = mycorr(xsol(cerebrum_rois, ind_lmci), abs(lmci_atrophy(cerebrum_rois)));
[Rall_ad, R_ad] = mycorr(xsol(cerebrum_rois,:), ad_atrophy(cerebrum_rois));
[~, ind_ad] = max(Rall_ad);
% [R_ad] = mycorr(xsol(cerebrum_rois, ind_ad), abs(ad_atrophy(cerebrum_rois)));

[Rall_korea_namci, R_korea_namci] = mycorr(xsol(ii,:), korea_namci_tau(ii));
[~, ind_korea_namci] = max(Rall_korea_namci);
% [R_korea_namci] = mycorr(xsol(ii, ind_korea_namci), abs(korea_namci_tau(ii)));
[Rall_korea_amci, R_korea_amci] = mycorr(xsol(ii,:), korea_amci_tau(ii));
[~, ind_korea_amci] = max(Rall_korea_amci);
% [R_korea_amci] = mycorr(xsol(ii, ind_korea_amci), abs(korea_amci_tau(ii)));
[Rall_korea_ad, R_korea_ad] = mycorr(xsol(ii,:), korea_ad_tau(ii));
[~, ind_korea_ad] = max(Rall_korea_ad);
% [R_korea_ad] = mycorr(xsol(ii, ind_korea_ad), abs(korea_ad_tau(ii)));

        figstr = [figuresavename '-emciatrophy'];
    yy =1*ballradius*emci_atrophy; %5
    [h, tiffname2] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_namcitau'];
    yy =ballradius*korea_namci_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname3] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-lmciatrophy'];
    yy =1*ballradius*lmci_atrophy; %5
    [h, tiffname4] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_amcitau'];
    yy =ballradius*korea_amci_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname5] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-adatrophy'];
    yy =1*ballradius*ad_atrophy; %5
    [h, tiffname6] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_adtau'];
    yy =ballradius*korea_ad_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname7] = glassbrain_display(yy(cerebrum_rois), figstr);

figure; 
cnt = 0;
for i = 1:length(tsamples)
    subplot(length(tsamples),5,cnt+1); subimage(imread(tiffnamex{i}));  title(['Model tau spread, t = ' num2str(tsamples(i))]); axis off
    cnt = cnt+5;
end
cnt = 0;
subplot(length(tsamples),5,cnt+2); h = plot(tsol, Rall_korea_namci, 'g-', tsol, Rall_korea_amci, 'r-', tsol, Rall_korea_ad, 'b-');  [h.LineWidth] = deal(2); 
title('Similarity index between model and tau SUVr'); xlabel('Model time (a.u.)'); ylabel('Similarity');
legend('naMCI', 'aMCI', 'AD');

subplot(length(tsamples),5,cnt+4); h = plot(tsol, Rall_emci, 'g-', tsol, Rall_lmci, 'r-', tsol, Rall_ad, 'b-');  [h.LineWidth] = deal(2); 
title('Similarity index between model and ADNI atrophy'); xlabel('Model time (a.u.)'); ylabel('Similarity');
legend('EMCI', 'LMCI', 'AD');

cnt = cnt + 5;
subplot(length(tsamples),5,cnt+2); subimage(imread(tiffname3));  title('naMCI tau SUVr'); axis off
subplot(length(tsamples),5,cnt+3); scatter(xsol(ii, ind_korea_namci), abs(korea_namci_tau(ii)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_korea_namci)]); xlabel('Model'); ylabel('naMCI tau SUVr');
subplot(length(tsamples),5,cnt+4); subimage(imread(tiffname2));  title('ADNI EMCI atrophy'); axis off
subplot(length(tsamples),5,cnt+5); scatter(xsol(cerebrum_rois, ind_emci), abs(emci_atrophy(cerebrum_rois)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_emci)]); xlabel('Model'); ylabel('ADNI EMCI atrophy');

cnt = cnt + 5;
subplot(length(tsamples),5,cnt+2); subimage(imread(tiffname5));  title('aMCI tau SUVr'); axis off
subplot(length(tsamples),5,cnt+3); scatter(xsol(ii, ind_korea_amci), abs(korea_amci_tau(ii)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_korea_amci)]); xlabel('Model'); ylabel('aMCI tau SUVr');
subplot(length(tsamples),5,cnt+4); subimage(imread(tiffname4));  title('ADNI LMCI atrophy'); axis off
subplot(length(tsamples),5,cnt+5); scatter(xsol(cerebrum_rois, ind_lmci), abs(lmci_atrophy(cerebrum_rois)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_lmci)]); xlabel('Model'); ylabel('ADNI LMCI atrophy');

cnt = cnt + 5;
subplot(length(tsamples),5,cnt+2); subimage(imread(tiffname7));  title('AD tau SUVr'); axis off
subplot(length(tsamples),5,cnt+3); scatter(xsol(ii, ind_korea_ad), abs(korea_ad_tau(ii)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_korea_ad)]); xlabel('Model'); ylabel('AD tau SUVr');
subplot(length(tsamples),5,cnt+4); subimage(imread(tiffname6));  title('ADNI AD atrophy'); axis off
subplot(length(tsamples),5,cnt+5); scatter(xsol(cerebrum_rois, ind_ad), abs(ad_atrophy(cerebrum_rois)), 20, 'b', 'fill'); lsline; 
title(['Peak R = ' num2str(R_ad)]); xlabel('Model'); ylabel('ADNI AD atrophy');


%% % Use later for plotting glassbrains of all tau groups 
%         figstr = [figuresavename '-patient_tau'];
%     yy =0.5*ballradius*patient_tau(:);
%     yy(Striatal_inds) = min(yy);
%     yy = min(yy, mean(yy) + 2*std(yy));
%     yy = max(0, yy - 1*std(yy));
%     [h, tiffname5] = glassbrain_display(yy, figstr);
%     [R3, p3] = corr(x0(ii), abs(patient_tau(ii)));
%     R3a = mycorr(xsol(ii,:), patient_tau(ii));
%     
%     figstr = [figuresavename '-koreapatient_tau1'];
%     yy =0.5*ballradius*korea_namci_tau;
%     yy(Striatal_inds) = min(yy);
%     yy = min(yy, mean(yy) + 2*std(yy));
%     yy = max(0, yy - 1*std(yy));
%     [h, tiffname6] = glassbrain_display(yy, figstr);
%     [R4, p4] = corr(x0(ii), abs(korea_namci_tau(ii)));
%     R4a = mycorr(xsol(ii,:), korea_namci_tau(ii));
% 
%     figstr = [figuresavename '-koreapatient_tau2'];
%     yy =0.5*ballradius*korea_amci_tau;
%     yy(Striatal_inds) = min(yy);
%     yy = min(yy, mean(yy) + 2*std(yy));
%     yy = max(0, yy - 1*std(yy));
%     [h, tiffname7] = glassbrain_display(yy, figstr);
%     [R5, p5] = corr(x0(ii), abs(korea_amci_tau(ii)));
%     R5a = mycorr(xsol(ii,:), korea_amci_tau(ii));
% 
%     figstr = [figuresavename '-koreapatient_tau3'];
%     yy =0.5*ballradius*korea_ad_tau;
%     yy(Striatal_inds) = min(yy);
%     yy = min(yy, mean(yy) + 2*std(yy));
%     yy = max(0, yy - 1*std(yy));
%     [h, tiffname8] = glassbrain_display(yy, figstr);
%     [R6, p6] = corr(x0(ii), abs(korea_ad_tau(ii)));
%     R6a = mycorr(xsol(ii,:), korea_ad_tau(ii));
% 
%     figure; 
%     subplot(5,2,1); subimage(imread(tiffname1)); title(' Early tau deposition sites, model prediction'); axis off
%     subplot(5,2,2); plot(tsol, R3a, 'g-', tsol, R4a, 'r:', tsol, R5a, 'r-', tsol, R6a, 'b-'); title('R, model versus patient AV1451-PET uptake (green: early patients, red: MCI, blue: AD');
%     subplot(5,2,3); subimage(imread(tiffname5)); title(' AV1451-PET uptake in early patients'); axis off
%     subplot(5,2,4); scatter(x0(ii), abs(patient_tau(ii)), 50, 'b', 'fill'); lsline; 
%     title(['Early tau model vs early patients AV1451-PET, R = ' num2str(R3)]); xlabel('Early tau model'); ylabel('Early patient group AV1451 uptake');
%     subplot(5,2,5); subimage(imread(tiffname6)); title(' AV1451-PET uptake in naMCI patients'); axis off
%     subplot(5,2,6); scatter(x0(ii), abs(korea_namci_tau(ii)), 50, 'b', 'fill'); lsline;
%     title(['Early tau model vs naMCI patients AV1451-PET, R = ' num2str(R4)]); xlabel('Early tau model'); ylabel('naMCI patient group AV1451 uptake');
%     subplot(5,2,7); subimage(imread(tiffname7)); title(' AV1451-PET uptake in aMCI patients'); axis off
%     subplot(5,2,8); scatter(x0(ii), abs(koreapatient_tau(ii,2)), 50, 'b', 'fill'); lsline;
%     title(['Early tau model vs aMCI patients AV1451-PET, R = ' num2str(R5)]); xlabel('Early tau model'); ylabel('aMCI patient group AV1451 uptake');
%     subplot(5,2,9); subimage(imread(tiffname8)); title(' AV1451-PET uptake in naMCI patients'); axis off
%     subplot(5,2,10); scatter(x0(ii), abs(koreapatient_tau(ii,3)), 50, 'b', 'fill'); lsline;
%     title(['Early tau model vs AD patients AV1451-PET, R = ' num2str(R6)]); xlabel('Early tau model'); ylabel('AD patient group AV1451 uptake');

    
    % % Stratify by model and amyloid positivity - failed, do not include
% R11 = mycorr(xsolnoint(cerebrum_rois,:), mciabneg_atrophy(cerebrum_rois));
% R12 = mycorr(xsolnoint(cerebrum_rois,:), mciabpos_atrophy(cerebrum_rois));
% R21 = mycorr(xsol(cerebrum_rois,:), mciabneg_atrophy(cerebrum_rois));
% R22 = mycorr(xsol(cerebrum_rois,:), mciabpos_atrophy(cerebrum_rois));
% 
% R11m = mycorr(xsolnoint(cerebrum_rois,:), mciabneg_met(cerebrum_rois));
% R12m = mycorr(xsolnoint(cerebrum_rois,:), mciabpos_met(cerebrum_rois));
% R21m = mycorr(xsol(cerebrum_rois,:), mciabneg_met(cerebrum_rois));
% R22m = mycorr(xsol(cerebrum_rois,:), mciabpos_met(cerebrum_rois));
% 
% R11ab = mycorr(fsolnoint(cerebrum_rois,:), mciabneg_abeta(cerebrum_rois));
% R12ab = mycorr(fsolnoint(cerebrum_rois,:), mciabpos_abeta(cerebrum_rois));
% R21ab = mycorr(fsol(cerebrum_rois,:), mciabneg_abeta(cerebrum_rois));
% R22ab = mycorr(fsol(cerebrum_rois,:), mciabpos_abeta(cerebrum_rois));
% 
% figure; 
% subplot(1,3,1); plot(tsolnoint, R11ab, 'r:', tsolnoint, R12ab, 'r-', tsol, R21ab, 'k:', tsol, R22ab, 'k-'); title('Empirical A\beta against ODE solution: the effect of amyloid-facilitation');
% legend('pure-tau model against amyloid-negative A\beta', 'pure-tau model against amyloid-positive A\beta', 'A\beta-facilitated tau model against amyloid-negative A\beta', 'A\beta-facilitated tau model against amyloid-positive A\beta');
% subplot(1,3,2); plot(tsolnoint, R11, 'r:', tsolnoint, R12, 'r-', tsol, R21, 'k:', tsol, R22, 'k-'); title('Empirical atrophy against ODE solution: the effect of amyloid-facilitation');
% legend('pure-tau model against amyloid-negative atrophy', 'pure-tau model against amyloid-positive atrophy', 'A\beta-facilitated tau model against amyloid-negative atrophy', 'A\beta-facilitated tau model against amyloid-positive atrophy');
% subplot(1,3,3); plot(tsolnoint, R11m, 'r:', tsolnoint, R12m, 'r-', tsol, R21m, 'k:', tsol, R22m, 'k-'); title('Empirical hypometabolism against ODE solution: the effect of amyloid-facilitation');
% legend('pure-tau model against amyloid-negative hypometabolism', 'pure-tau model against amyloid-positive hypometabolism', 'A\beta-facilitated tau model against amyloid-negative hypometabolism', 'A\beta-facilitated tau model against amyloid-positive hypometabolism');

%% ***************************************************************************************
% Optimization, repeated seeding, network scrambling, etc
% This section runs above code and functions multi9ple times, wherther
% repeated seedings et cor within an optimization procedure
% Only uncomment these sections below if needed

%% ADDED NEW 11/20/19: repeat over all possible seed locations (single seeds)

% % Doing this for aMCI subjects only - the most clinically interesting group
% % Keep this commneted out for routine runs - this takes FOREVER
% 
%     xdrive_orig = xdrive;
%     seedvec = [1:34, 69:77; 35:68, 78:86];
%     fdrive = 0.5*nc_met.*nc_app;
%     f0 = zeros([size(fdrive), M]);
%     x0 = zeros([size(xdrive), M]);
%     xf0 = [ x0(:); f0(:)];  
%     for q1 = 1:size(seedvec, 2)      
%         xdrive = zeros(nroi,1);
%         xdrive(seedvec(:,q1)) = 800;
% 
%         [tsol, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%         nt = size(xfsol_all,1);
%         xsol_all = xfsol_all(:, 1:nroi*M); 
%         xsol_all = reshape(xsol_all, [nt, nroi,M]);
%         xsol_all = permute(xsol_all, [2,1,3]);  
%         xsol = xsol_all(:,:,M);
%         tmp = mycorr(xsol(ii,:), korea_amci_tau(ii)); 
%         RKorea_amci_repeatseeding(q1) = max(tmp);
%     end
%     xdrive = xdrive_orig;
% save RKorea_amci_repeatseeding.mat RKorea_amci_repeatseeding;
% figure; bar(RKorea_amci_repeatseeding); title('Rmax of each region being seeded in turn');


%% ADDED NEW 3/20/19: explore model behaviour against production rate, etc
% %     alphavec = linspace(0.01, 0.1, 10);
% origalpha = alpha;
% alphavec = logspace(-3, -1, 12);
% trange = [0, 200];
% tsol = cell(length(alphavec),1);
% total_tau = cell(length(alphavec),1);
% RKorea_amci_all = cell(length(alphavec),1);
% total_tau_max = []; 
% RKorea_amci = [];
% 
% for q1 = 1:length(alphavec)
%     alpha = alphavec(q1); 
%     [tsol{q1}, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%     nt = size(xfsol_all,1);
%     xsol_all = xfsol_all(:, 1:nroi*M); 
%     xsol_all = reshape(xsol_all, [nt, nroi,M]);
%     xsol_all = permute(xsol_all, [2,1,3]);  
%     xsol = xsol_all(:,:,M);
%     total_tau{q1} = sum(xsol,1); 
%     total_tau_max(q1) = max(total_tau{q1});
%     tmp = mycorr(xsol(ii,:), korea_amci_tau(ii)); 
%     RKorea_amci_all{q1} = tmp;
%     RKorea_amci(q1) = max(tmp);
% end
% alpha = origalpha;
% figure; subplot(221); hold on;
% for q1 = 1:length(alphavec)
%     plot(tsol{q1}, total_tau{q1},'LineWidth', 3); 
% end
% title('Brain-wide tau burden for varying monomer production rates'); xlabel('Model time (a.u.)'); ylabel('AND predicted tau burden (a.u.)');
% subplot(222); plot(alphavec, total_tau_max,'LineWidth', 3);  title('Tau burden vs. monomer production rate'); xlabel('monomer production rates (a.u.)'); ylabel('AND predicted tau burden (a.u.)'); 
% 
% subplot(223); hold on;
% for q1 = 1:length(alphavec)
%     plot(tsol{q1}, RKorea_amci_all{q1},'LineWidth', 3); 
% end
% title(sprintf('Correaltion of regional tau SUVR and AND model at different monomer production rates')); xlabel('Model time (a.u.)'); ylabel(' R');
% subplot(224); plot(alphavec, RKorea_amci,'LineWidth', 3);  title('Peak correlation vs. monomer production rate'); xlabel('monomer production rate (a.u.)'); ylabel('Peak R'); 
% % cftool(alphavec, total_tau_max);
% 
% % combo plot for varying diffusivities
% figure; subplot(121); hold on; title('Tau burden vs. monomer production rate'); xlabel('monomer production rates (a.u.)'); ylabel('AND predicted tau burden (a.u.)'); 
% subplot(122); hold on; title('Model time of peak Tau burden vs. monomer production rate'); xlabel('monomer production rates (a.u.)'); ylabel('Model time of peak tau burden (a.u.)'); 
% 
% origbeta = beta;
% betavec = origbeta*logspace(-1, 1, 3);
% 
% for q2 = 1:length(betavec)
%     beta = betavec(q2); 
%     tsol = cell(length(alphavec),1);
%     total_tau = cell(length(alphavec),1);
%     RKorea_amci_all = cell(length(alphavec),1);
%     for q1 = 1:length(alphavec)
%         alpha = alphavec(q1);  % EC seeding    
%         [tsol{q1}, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%         nt = size(xfsol_all,1);
%         xsol_all = xfsol_all(:, 1:nroi*M); 
%         xsol_all = reshape(xsol_all, [nt, nroi,M]);
%         xsol_all = permute(xsol_all, [2,1,3]);  
%         xsol = xsol_all(:,:,M);
%         total_tau{q1} = sum(xsol,1);
%         total_tau_max(q1) = max(total_tau{q1});
%         ind = find(total_tau{q1}>0.99*total_tau_max(q1), 1);
%         tmx(q1) = tsol{q1}(ind);
%     end
%     subplot(121); plot(alphavec, total_tau_max,'LineWidth', 3);
%     subplot(122); plot(alphavec, tmx,'LineWidth', 3);
% end
% alpha = origalpha;
% beta = origbeta;

%% ADDED NEW 5/24/18: optimize model parameters

% Keep this commneted out for routine runs - this takes FOREVER
% 
%     sM_orig = sM;
%     ctau_orig = ctau;
%     alpha_orig = alpha;
%     beta = 20;
%     gama = 0; % enforce no interaction
%     lambda = 15;    
%     sMvec = [0.5, 1,2,3]; %[1,2];
%     ctauvec = [2, 5, 10, 15];%[10,15]; 
%     alphavec = [0.01, 0.02, 0.05, 0.25]*beta;
%     
%     for q1 = 1:length(sMvec)
%         sM = sMvec(q1);
%         %     A = ones(M,M)/sM;  % coeffs a_ij, change later
% %             A = sM*sM./((1:M).'*(1:M));  % coeffs a_ij= 1/ij, based on statistical mechanics
%               A = ((1:M).'*(1:M))/sM/sM .* exp(-(1:M).'*(1:M)/sM/sM);  % coeffs a_ij= 1/ij, based on statistical mechanics
%         %     dconst = sM/M;   % parameter of exponential decay mode governing the diffusion constants of oligomers of bvarying length
% %             diffvec = exp(-(1:M).'/sM);
% %             diffvec =  2*sM*sM*((1:M).^(-2)).';
%               diffvec = (1:M).'/sM .* exp(-(1:M).'/sM);
%         for q2 = 1:length(ctauvec)
%             ctau = ctauvec(q2);
%             for q3 = 1:length(alphavec)
%                 alpha = alphavec(q3);
%                 [tsol, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%                 nt = size(xfsol_all,1);
%                 xsol_all = xfsol_all(:, 1:nroi*M); 
%                 xsol_all = reshape(xsol_all, [nt, nroi,M]);
%                 xsol_all = permute(xsol_all, [2,1,3]);  
%                 xsol = xsol_all(:,:,M);
% 
%                 tmp = mycorr(xsol(ii,:), korea_amci_tau(ii)); 
%                 RKorea_amci(q1,q2,q3) = max(tmp);
%             end
%         end
%     end
%     RKorea_amci,
%     [Rmax, ind] = max(RKorea_amci(:));
%     [q1, q2, q3] = ind2sub([length(sMvec), length(ctauvec), length(alphavec)], ind);
%     Rmax,
%     sM_opt = sMvec(q1),
%     ctau_opt = ctauvec(q2),
%     alpha_opt = alphavec(q3),
%     save RKorea_amci.mat RKorea_amci;
%     figure; surf(alphavec, sMvec, squeeze(RKorea_amci(:,q2,:)));
%     sM = sM_orig;
%     ctau = ctau_orig;
%     alpha = alpha_orig;



%% 8) Node scrambling

% 8a) Node scrambling for tau

% if edge_scrambling
%     ntrials = 2000; 
%     Rqsave = zeros(ntrials,1);
%     for i = 1:ntrials        
%         Crand = C(:);
%         randinds = randperm(length(Crand));
%         Crand = Crand(randinds);
%         Crand = reshape(Crand, nroi, nroi);
%         Crand = tril(Crand,-1) + triu(Crand.',1);
%         rowdegree = (sum(Crand, 2)).';
%         coldegree = sum(Crand, 1);
%         Lret = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Crand; % * diag(1./(sqrt(coldegree)+eps)) ;
%         [tsol, xfsol_all] = ode45(@rajfranchi_ode, trange, xf0, odeopts);
%         nt = size(xfsol_all,1);
%         xsol_all = xfsol_all(:, 1:nroi*M); 
%         xsol_all = reshape(xsol_all, [nt, nroi,M]);
%         xsol_all = permute(xsol_all, [2,1,3]);  
%         xsol = xsol_all(:,:,M);
% %         total_tau{i} = sum(xsol,1); 
% %         total_tau_max(i) = max(total_tau{i});
%         tmp = mycorr(xsol(ii,:), korea_amci_tau(ii)); 
%         Rqsave(i) = max(tmp);
%     end
% 
% %         plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-', tsol, R21, 'b:', tsol, R22, 'b-'); 
% %     title('R, model versus empirical AD atrophy and its randomized node permutations'); hold on;
% %     legend( 'EMCI', 'LMCI', 'AD', 'MCI- A\beta-', 'MCI- A\beta+', 'permuted AD');
%     figure;
%     hist(Rqsave, 30, 'FaceColor', 'g'); hold on; 
%      title('Peak R of randomized connectome permutations, evaluated on the Korea aMCI cohort');
%     hold off;
%      
%     %  restore correct C etc
%     Cret = C; %retro
%     rowdegree = (sum(Cret, 2)).';
%     coldegree = sum(Cret, 1);
%     Lret = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Cret; % * diag(1./(sqrt(coldegree)+eps)) ;
% 
% end

% 8b) Node scrambling
% if node_scrambling
%     ntrials = 1000; 
%     Rqsave = zeros(ntrials,1);
%     q = ad_atrophy(cerebrum_rois);
%     figure;
%     subplot(1,2,1); plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-', tsol, R21, 'b:', tsol, R22, 'b-'); 
%     title('R, model versus empirical AD atrophy and its randomized node permutations'); hold on;
%     for i = 1:ntrials
%         q = q(randperm(length(cerebrum_rois)));
%         Rq = mycorr(xsol(cerebrum_rois,:), q);
%     %     Rq = overlap_measure(xsol(cerebrum_rois,:), q, pthr); 
%         if mod(i,5)==0, plot(tsol, Rq, 'g-'); end
%         Rqsave(i) = max(Rq(:));
%         
%     end
%     legend( 'EMCI', 'LMCI', 'AD', 'MCI- A\beta-', 'MCI- A\beta+', 'permuted AD');
%     subplot(1,2,2);
%     hist(Rqsave, 20, 'FaceColor', 'g'); hold on; 
%     line([max(Rad), max(Rad)], [0, ntrials/10], 'Color', 'k'); 
%     line([max(Remci), max(Remci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', ':'); 
%     line([max(Rlmci), max(Rlmci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', '-'); 
%     line([max(R21), max(R21)], [0, ntrials/10], 'Color', 'b', 'LineStyle', ':'); 
%     line([max(R22), max(R22)], [0, ntrials/10], 'Color', 'b', 'LineStyle', '-'); 
%      title('Peak R of empirical regional atrophy and randomized permutations');
%      legend('permuted AD', 'AD', 'EMCI', 'LMCI', 'MCI- A\beta-', 'MCI- A\beta+');
% 
%     % Do same for Abeta dist
%     Remci = mycorr(fsol(cerebrum_rois,:), emci_abeta(cerebrum_rois));
%     Rlmci = mycorr(fsol(cerebrum_rois,:), lmci_abeta(cerebrum_rois));
%     Rad = mycorr(fsol(cerebrum_rois,:), ad_abeta(cerebrum_rois));
%     Rqsave = zeros(ntrials,1);
%     q = ad_abeta(cerebrum_rois);
%     figure;
%     subplot(1,2,1); plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-', tsol, R21ab, 'b:', tsol, R22ab, 'b-'); 
%     title('R, model versus empirical AD A\beta distribution and its randomized node permutations'); hold on;
%     for i = 1:ntrials
%         q = q(randperm(length(cerebrum_rois)));
%         Rq = mycorr(fsol(cerebrum_rois,:), q); 
%     %     Rq = overlap_measure(fsol(cerebrum_rois,:), q, pthr); 
%         if mod(i,5)==0, plot(tsol, Rq, 'g-'); end
%         Rqsave(i) = max(Rq(:));
%     end
%     legend( 'EMCI', 'LMCI', 'AD', 'MCI- A\beta-', 'MCI- A\beta+', 'permuted AD');
%     subplot(1,2,2);
%     hist(Rqsave, 20, 'FaceColor', 'g'); hold on; 
%     line([max(Rad), max(Rad)], [0, ntrials/10], 'Color', 'k'); 
%     line([max(Remci), max(Remci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', ':'); 
%     line([max(Rlmci), max(Rlmci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', '-'); 
%     line([max(R21ab), max(R21ab)], [0, ntrials/10], 'Color', 'b', 'LineStyle', ':'); 
%     line([max(R22ab), max(R22ab)], [0, ntrials/10], 'Color', 'b', 'LineStyle', '-'); 
%      title('Peak R of empirical regional A\beta and randomized permutations');
%      legend('permuted AD', 'AD', 'EMCI', 'LMCI', 'MCI- A\beta-', 'MCI- A\beta+');
% end
% 
%  % Plots of global pathology burden
%  for i = 1:length(tsol)
%     q = xsol(:,i);
%     globaltau(i) = sum(q);
%     temptau(i) = sum(q(lobes==4));
%     q = fsol(:,i);
%     globalab(i) = sum(q);
%     tempab(i) = sum(q(lobes==4));
%  end
%  for i = 1:length(tsolnoint)
%     q = xsolnoint(:,i);
%     globaltaunoint(i) = sum(q);
%     temptaunoint(i) = sum(q(lobes==4));
%  end
%  figure; plot(tsol, globalab, 'r-', tsol, tempab, 'r:', tsol, globaltau, 'b-', tsol, temptau, 'b:', tsolnoint, globaltaunoint, 'g-', tsolnoint, temptaunoint, 'g:');
%  title('Global and temporal burden of modeled A\beta and tau');


%% % Glassbrains of ODE dyanmics
% figure;      
% for i = 1:length(tsamples)
%     figstr = [figuresavename '-ODEx' num2str(i)];
%     [dummy, ind] = min(abs(tsol-tsamples(i)));
%     yy = 10*ballradius*xsol(:,ind);
%     yy = min(yy, mean(yy) + 1*std(yy));
%     %yy = max(0, yy - 1*std(yy));
%     [h, tiffnamex{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
%     figstr = [figuresavename '-ODEf' num2str(i+length(tsamples))];
%     yy = 10*ballradius*fsol(:,ind);
%     yy = min(yy, mean(yy) + 1*std(yy));
%     %yy = max(0, yy - 1*std(yy));
%     [h, tiffnamef{i+length(tsamples)}] = glassbrain_display(yy(cerebrum_rois), figstr);
%     
%     figstr = [figuresavename '-ODExnoint' num2str(i)];
%     [dummy, ind] = min(abs(tsolnoint-tsamples(i)));
%     yy = 10*ballradius*xsolnoint(:,ind);
%     yy = min(yy, mean(yy) + 1*std(yy));
%     %yy = max(0, yy - 1*std(yy));
%     [h, tiffnamexnoint{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
% %     figstr = [figuresavename '-ODEfnoint' num2str(i+length(tsamples))];
% %     yy = 10*ballradius*fsolnoint(:,ind);
% %     yy = min(yy, mean(yy) + 1*std(yy));
% %     %yy = max(0, yy - 1*std(yy));
% %     [h, tiffnamefnoint{i+length(tsamples)}] = glassbrain_display(yy(cerebrum_rois), figstr);
% end
% figure; 
% cnt = 0;
% for i = 1:length(tsamples)
%     subplot(length(tsamples),3,cnt+1); subimage(imread(tiffnamef{i+length(tsamples)}));  title(['A\beta spread, t = ' num2str(tsamples(i))]); axis off
%     subplot(length(tsamples),3,cnt+2); subimage(imread(tiffnamex{i}));  title(['tau spread, t = ' num2str(tsamples(i))]); axis off
%     subplot(length(tsamples),3,cnt+3); subimage(imread(tiffnamexnoint{i}));  title(['tau spread, t = ' num2str(tsamples(i))]); axis off
%     cnt = cnt+3;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%  %%%%%%%%   INTERNAL FUNCTIONS  %%%%%%%%%%%

    function xfprime = myode(t, xf)
        beta = 2; %1 .5
        alpha = 2; %*beta; %1.5 3 1
        gama = 0.5; %*beta; %0.5
        lambda = 35; % 20 15 10
        eta = 4;
        xthr = alpha*mean(x0)*tmax/8;
        fthr = alpha*mean(f0)*tmax/8;
        
        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
        gfun = (t/lambda)*exp(-t/lambda);
%         gfun = 5/tmax;    % 5 for tau, 20 for amyloid
%         xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0;
        xfprime(1:nroi) =        -beta*Lret*xu  + alpha*gfun*x0 + eta*max(xu-xthr, 0); % + gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0 + eta*max(fu-fthr, 0) ;
        xfprime = xfprime(:);
    end

    function xfprime = myode_nointeraction(t, xf)
        beta = 1; %1 .5
        alpha = 1; %*beta;  %1.5 3 1
        lambda = 15; % 20 10
        eta = 0.75;
        xthr = alpha*mean(x0)*tmax/5;
        fthr = alpha*mean(f0)*tmax/5;

        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
        gfun = (t/lambda)*exp(-t/lambda);
%         xfprime(1:nroi) =        -beta*Lret*xu  + alpha*gfun*x0; % + gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0;
        xfprime(1:nroi) =          -beta*Lret*xu  + alpha*gfun*x0 + eta*max(xu-xthr, 0) ; %+ gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi)   = -beta* Lret*fu + alpha*gfun*f0 + eta*max(fu-fthr, 0) ;
        xfprime = xfprime(:);
    end

    function xfprime = myode2(t, xf)
        beta = 4;
        alpha = 0.1*beta;
        gama = 0.3*beta;
        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
% %         xfprime(1:nroi) =      -beta*Lret*xu  +alpha*(fu.*xu);
%         xfprime(1:nroi) =        -beta*Lret*xu  + alpha*gfun*x0 + gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0;
        xfprime(1:nroi) =          -beta*Lret*xu  + alpha* tau_production_fun(xu) + gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi) =   -beta* Lret*fu + alpha* Abeta_production_fun(fu);
        xfprime = xfprime(:);
        
        function hfun = tau_production_fun(x)
            sigm = max(x0);  %0.08
            %hfun = erf((x-sigm)/sigm);
            hfun = 0.01*sqrt(abs(x)/sigm);
        end
        function hfun = Abeta_production_fun(f)
            sigm = max(f0);  %0.03;
            %hfun = erf((f-sigm)/sigm);
            hfun = 0.01*sqrt(abs(f)/sigm);
        end
    end

    function [R, Rmax] = mycorr(xx, yy)
        if ~exist('corrtype', 'var') || strcmp(corrtype, 'Pearson') 
            R = corr(xx,yy);
        elseif strcmp(corrtype, 'Spearman')
            inds = (yy>median(yy(:)));
            xx = xx(inds,:);
            yy = yy(inds);
            R = corr(xx,yy, 'type', 'Spearman');
        elseif strcmp(corrtype, 'Dotproduct')
            nvecs = size(xx,2);
            for kk = 1:nvecs
                xx1 = xx(:,kk)/norm(xx(:,kk));
                R(kk) = xx1.'*yy;
            end
        elseif strcmp(corrtype, 'Overlap')
            % new measure of correlation between model and empirical
            % assumes that only xx can have multipel columns
            nvecs = size(xx,2);
            for kk = 1:nvecs
                x = xx(:,kk);
                y = yy;
                if pthr==0.5
                    thrx= median(x);
                    thry= median(y);
                else
                    sx = sort(x, 'Ascend');
                    thrx = sx(round(pthr*length(sx)));
                    sy = sort(y, 'Ascend');
                    thry = sy(round(pthr*length(sy)));
                end
                x1 = double(x > thrx);
                y1 = double(y > thry);
                R(kk) = sum(x1.*y1)/length(x);
            end
        elseif strcmp(corrtype, 'boxcox')
            % use boxcox to resample to match distributions
            bcoxfact = 1; %0.5;
            yy1 = boxcox(yy(:));
            [transdat, lambda] = boxcox(xx(:)) ; 
            lambda1 = bcoxfact*lambda;
            x = boxcox(lambda1, xx(:));
            xx1 = reshape(x, size(xx));
            R = corr(xx1, yy1);
        end
        R = R(:);
        Rmax = max(R);
    end
    
    function xfprime = rajfranchi_ode(t, xf)
        xu = xf(1:nroi*M);
        fu = xf(nroi*M+1:2*nroi*M);
        xu = reshape(xu, [nroi,M]); 
        fu = reshape(fu, [nroi,M]); 
        xprime = zeros(size(xu));
        fprime = zeros(size(fu));            
        gfun = (t/lambda)*exp(-t/lambda);

        for m = 1:M-1
            ss = cAb* fu*(A(m,:)).'; 
            fprime(:,m) = -beta*diffvec(m)*Lret*fu(:,m) - fu(:,m).*ss ;
        end
        fprime(:,1) = fprime(:,1)  + alpha*gfun*fdrive;
        %fprime(:,1) = fprime(:,1)  + alpha*gfun*fu(:,1) ;
        fprime(:,M) = -beta*diffvec(M)*Lret*fu(:,M) ;  % at M you have plaque/tangle, hence it never decreases

        for m = 2:M
            q = 0;
            for k = 1:m-1
                q = q+0.5*cAb*A(k,m-k)*fu(:,k).*fu(:,m-k);
            end
            fprime(:,m) = fprime(:,m) + q;
        end

        for m = 1:M-1
            ss = ctau* xu*(A(m,:)).'; 
            xprime(:,m) = -beta*diffvec(m)*Lret*xu(:,m) - xu(:,m).*ss ;
        end
        q = fu*diffvec;
        xprime(:,1) = xprime(:,1)  + alpha*gfun*xdrive + gama*(q.*xu(:,1)) ;
        %xprime(:,1) = xprime(:,1)  + alpha*gfun*xu(:,1) + gama*(q.*xu(:,1)) ;
        xprime(:,M) = -beta*diffvec(M)*Lret*xu(:,M) ;  % at M you have plaque/tangle, hence it never decreases
        for m = 2:M
            q = 0;
            for k = 1:m-1
                q = q+0.5*ctau*A(k,m-k)*xu(:,k).*xu(:,m-k);
            end
            xprime(:,m) = xprime(:,m) + q;
        end
        xfprime = [xprime, fprime];
        xfprime = xfprime(:);
    end

function [ha, tiffname] = glassbrain_display(yy, figstr)

% Glass brain Display and figure save settings
    customLUT='hot'; % '' or e.g. 'winter', 'jet';
    pipesize = 'Normal'; % 'Small','Large','XL'
    nodesize = 'XL';
    atlas_type = 'MS86'; %'MS86';  % or 'ADNI116'
    baseStructFileName = ['brainography-master-Ashish' filesep 'sample_files' filesep 'ashish_86_region_base.mat']; %Need to be in local path or set absolute path in this string
    gbplottype = 'glassbrain'; % 'glassbrain'; % or 'surface'
    pipesonoff = 1; % or 1
    connectivityMatrix = [];
    
    colstr = {'c-', 'r-', 'g-', 'bl-', 'b:'};
    map = colormap('jet');
    fig_fontsize = 12;  % font within matlab figures
    fig_font = 'Cambria';
    fig_linewidth = 2;

    % main plotting
    v_all = yy(:);
    lq = mean(v_all) - 1*std(v_all); uq = mean(v_all) + 2*std(v_all);
    yy = max(yy, lq); yy = min(yy, uq); 
    yy = yy - min(yy(:)); 
    %yy = yy/max(yy(:));
    
    tiffname = [figstr '.tif'];
    v = yy;
    regionvalues = v(cerebrum_rois);
    figure;
    
    if suppress_glassbrain
        ha = gcf;
        export_fig(tiffname,'-tif','-r300','-a2');
        close(ha);
        return; 
    end

    H = brainomatic(baseStructFileName,atlas_type, gbplottype,figstr,regionvalues,connectivityMatrix,pipesonoff,customLUT,pipesize,nodesize);
    figstr = {[figstr '_axial.png'], [figstr '_sagittal.png']};
    imdisp(figstr, 'Size', [1,2], 'Border', [0, 0]);
    %     figure; montage(figstr,  'Size', [1 NaN]);  %axis normal;
    ha = gcf;
%     set(ha, 'Units', 'normalized', 'OuterPosition', [0, 0, 0.8, 0.3], 'Position', [0, 0, 0.8, 0.3]);
%     set(ha, 'Units', 'pixels');
    title(sprintf('Graph eigenmode'), 'FontSize', 16, 'FontName', 'Cambria');
%     ha = tightfig(ha);
    set(ha, 'Color', 'w');
    export_fig(tiffname,'-tif','-r300','-a2');
    close(ha);
end

end