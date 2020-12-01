function TauAmyloid_Joint_eNDM_Git()

% Implements the joint tau amyloid co-evolving eNDM model, as presented in the paper Raj et al, Nature Communications 2020
% x is the modeled tau, f is modeled amyloid. 
% production term or drive for x is seeded at entorhinal cortex (EC) unless otherwise required. In that case change seeds variable accordingly. 
% Uses dynamics driven by x0, f0 (production terms)
% model is evaluated and tested on empirical regional data from ADNI. In
% this code only group average regional data are used. Consequently
% individualized fitting of model parameters was not attempted or included
% here. Brute fornce grid search was used on group data to optimize model
% parameters - but this proceudre is not included in this code. If
% inetersted in re-optimizing, run the internal function that generates the
% model (see below) within nested for loops.

% 
% Revision of NComms paper: Added several new simulations in reponse to
% NComms reviews. These include new spread modes (Conn, DZistance, etc) and
% new interaction models (no-int, 1-way int, 2-way int)
% Also added new statistical readouts (AIC, Fisher)

use_conn = 'Conn'; % 'Conn'; % 'FiberDistance'  'EuclideanDistance' 'ComboCandD'
% Note: ComboCandD is yet to be properly explored; Conn is better than other models

% Select which connectome to use (here parameters were optimized over Cornell connectome). If using HCP, may need to re-optimize parameters
    use_hcp = 0;
    use_Cornell = 1;
% parameters for regional data rescaling
    sig = 2;
    sig_abeta = 2;
    sig_tau = 2;
    rescale_method = 'unit'; %'logistic';
    
% variables for glassbrain plotting
    figuresavename = 'test1';
    suppress_glassbrain = 1; % use this to suppress the (very slow) glassbrain rendering process
    ballradius = 2;  % just for glassbrain plotting
    fontsz = 16; % contrl sfont size on figs and plots

% Variables that govern ODE solver
    tmax = 30; % maximum time range for simulations
    trange = [0,tmax];
    tsamples = [5 10 15 20];  % sampled time pints for glass brain plotting

% permutation testing of ODE
    node_scrambling = 0;  
    edge_scrambling = 0; 
    ntrials = 1000; 
    repeat_seeding = 0;
corrtype = 'Pearson'; % 'Pearson' 'Spearman' or 'Overlap' or 'Dotproduct' or 'boxcox'; was DotProduct
pthr = 0.5; % only used for corrtype overlap; higher is more AD specific

% load and organize connectome data
Ddk = csvread('mean80_fiberlength.csv', 0, 0);     
S = load('Euclidean_dist_matrix.mat');
Ddk_euc = S.Euclidean_dist_matrix;
s = load ('FSatlas_lobes.mat');
lobes = s.lobes;
if use_hcp
    permHCP = [19:52, 53:86, 1:9, 10:18];  % permutes 86 x 86 HCP conn matrix so that subcorts are at the end, as in previous work
    Cdk_hcp = csvread('mean80_fibercount.csv', 1, 0);         
    Cdk_hcp = Cdk_hcp(permHCP, permHCP);
    C = Cdk_hcp;
elseif use_Cornell       
    permEve = [1:68, 78:86, 69:77];  % permutes 86 x 86 Cornell conn matrix so that subcorts are at the end, as in previous work
    load meanACS69.mat meanACS;
    Cdk_eve = meanACS;
    Cdk_eve = Cdk_eve(permEve, permEve);
    C = Cdk_eve;
end
nroi = size(C,1);

switch use_conn 
    case 'Conn'
       C=C;
    case 'FiberDistance'      
        ss = std(Ddk(:));
        Ddk(Ddk==0) = 1e6;
    %     C = exp(-Ddk/ss);
        C = exp(-Ddk.*Ddk/2/ss/ss);
        C = C.*(1-diag(ones(nroi,1)));
    case 'ComboCandD'
        ss = std(Ddk(:));
        Ddk(Ddk==0) = 1e6;
    %     C1 = exp(-Ddk/ss);
        C1 = exp(-Ddk.*Ddk/2/ss/ss);
        C1 = C1.*(1-diag(ones(nroi,1)));
        C = C.*C1;
    case 'EuclideanDistance'
        ss = std(Ddk_euc(:));
    %     C = exp(-Ddk_euc/ss);
        C = exp(-Ddk_euc.*Ddk_euc/2/ss/ss);
        C = C.*(1-diag(ones(nroi,1)));
end


% Define some useful indices
ECLind = 5; ECRind = 34+ECLind;  
ITLind = 8; ITRind = 34+ITLind;  
AMLind = 75; AMRind = 84; 
HPLind = 74; HPRind = 83; 
CDLind = 71; CDRind = 80; 
PULind = 72; PURind = 81;
latOCCLind = 10; latOCCRind = 44;
ThLind = 70; ThRind = 79;
Striatal_inds = [70 71 72 73 76 79 80 81 82 85]; % included thalamus
% Striatal_inds = [71 72 73 76 80 81 82 85];  
cerebrum_rois = 1:86;
noncerebrum_rois = 87:90;
cerebellum_inds = [69,78];
ii = setdiff(cerebrum_rois, Striatal_inds);  % all regions not in striatum or cerebellum

seeds = [ECLind, ECRind];  % default
% seeds = [HPLind, HPRind];  
% seeds = [ITLind, ITRind];  


nroi = size(C,1);

left_rois = [1:34, 69:77];
right_rois = [35:68, 78:86];

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


%% Load data

% APP map, artificially created
% (occ and parital are 0.5 of others, striatals are 1.1x)

nc_app = ones(nroi, 1);
nc_app(lobes==3) = 0.5;
nc_app(lobes==1) = 1.1;
nc_app([CDLind, CDRind, PULind, PURind]) = 1.1;
% MAPT map, artificially created (occ and striatal are 0.5 of others)
nc_mapt = ones(nroi, 1);
nc_mapt(lobes==3) = 0.5;
nc_mapt([CDLind, CDRind, PULind, PURind]) = 0.5;

% Rest of ADNI group data are contained in the .mat file
load ADNI_group_data.mat nc_met emci_met lmci_met ad_met nc_abeta emci_abeta lmci_abeta ad_abeta nc_tau emci_tau lmci_tau ad_tau emci_atrophy lmci_atrophy ad_atrophy;

% rescale brain data
temp = rescale_braindata([nc_met emci_met lmci_met ad_met], sig, rescale_method, 'allcols');
nc_met = temp(:,1);
emci_met = temp(:,2);
lmci_met = temp(:,3);
ad_met = temp(:,4);
% Note 5 factor below......
temp = 5*rescale_braindata([nc_abeta emci_abeta lmci_abeta ad_abeta], sig_abeta, rescale_method, 'allcols');
nc_abeta = temp(:,1);
emci_abeta = temp(:,2);
lmci_abeta = temp(:,3);
ad_abeta = temp(:,4);
temp = rescale_braindata([nc_tau emci_tau lmci_tau ad_tau], sig_tau, rescale_method, 'allcols');
nc_tau = temp(:,1);
emci_tau = temp(:,2);
lmci_tau = temp(:,3);
ad_tau = temp(:,4);

temp = rescale_braindata([emci_atrophy lmci_atrophy ad_atrophy], sig, rescale_method, 'allcols');
emci_atrophy = temp(:,1);
lmci_atrophy = temp(:,2);
ad_atrophy = temp(:,3);


% DEBGUG
figure; 
subplot(3,1,1); bar(lmci_tau); title('LMCI Tau');
subplot(3,1,2); bar(lmci_atrophy); title('LMCI atrophy');
subplot(3,1,3); bar(lmci_abeta); title('LMCI A-beta');

% Plot some relevant empirical regional patterns
    figstr = [figuresavename '-ncmet'];
        yy =2*ballradius*nc_met;
    [h, tiffname1] = glassbrain_display(yy(cerebrum_rois), figstr);
    figstr = [figuresavename '-ncapp'];
        yy =2*ballradius*nc_app;
    [h, tiffname2] = glassbrain_display(yy(cerebrum_rois), figstr);
    figstr = [figuresavename '-ncmapt'];
        yy =2*ballradius*nc_mapt;
    [h, tiffname3] = glassbrain_display(yy(cerebrum_rois), figstr);
    figure; 
    subplot(1,3,1); subimage(imread(tiffname1));  title(' Healthy baseline metabolic load '); axis off
    subplot(1,3,2); subimage(imread(tiffname2));  title(' Healthy APP gene expression '); axis off
    subplot(1,3,3); subimage(imread(tiffname3));  title(' Healthy MAPT gene expression '); axis off

%% % 6) Solve the differential equations numerically

% x = tau, f = Abeta

x0 = zeros(nroi,1);
x0(seeds) = 1;
x0 = x0/norm(x0); 
f0 = nc_met.*nc_app; 
f0 = f0/norm(f0);
xf0 = [zeros(size(nc_met)); zeros(size(nc_met))];  
odeopts = odeset('NonNegative', 1:length(xf0));
[tsol, xfsol] = ode45(@myode, trange, xf0, odeopts);
% [tsol, xfsol] = ode45(@myode2, trange, xf0, odeopts); % 2-way interaction model
xsol = xfsol(:, 1:nroi); xsol = xsol.';
fsol = xfsol(:, nroi+1:2*nroi); fsol = fsol.';
[tsolnoint, xfsolnoint] = ode45(@myode_nointeraction, trange, xf0, odeopts);  % no interaction model
xsolnoint = xfsolnoint(:, 1:nroi); xsolnoint = xsolnoint.';
fsolnoint = xfsolnoint(:, nroi+1:2*nroi); fsolnoint = fsolnoint.';


%% % 6-7a) Plot Rt curves w.r.t. ADNI regional data

% First the ODE dynamics of fsol and empirical Abeta
    
for i = 1:length(tsamples)
    figstr = [figuresavename '-ODEf' num2str(i)];
    [dummy, ind] = min(abs(tsol-tsamples(i)));
    yy = 4*ballradius*fsol(:,ind);  %5
%     yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamef{i}] = glassbrain_display(yy(cerebrum_rois), figstr); 
end
[Rall_emci, R_emci] = mycorr(fsol(cerebrum_rois,:), emci_abeta(cerebrum_rois));
% [~, ind_emci] = max(Rall_emci);
% [R_emci] = mycorr(fsol(cerebrum_rois, ind_emci), abs(emci_abeta(cerebrum_rois)));
[Rall_lmci, R_lmci] = mycorr(fsol(cerebrum_rois,:), lmci_abeta(cerebrum_rois));
% [~, ind_lmci] = max(Rall_lmci);
% [R_lmci] = mycorr(fsol(cerebrum_rois, ind_lmci), abs(lmci_abeta(cerebrum_rois)));
[Rall_ad, R_ad] = mycorr(fsol(cerebrum_rois,:), ad_abeta(cerebrum_rois));
% [~, ind_ad] = max(Rall_ad);
% [R_ad] = mycorr(fsol(cerebrum_rois, ind_ad), abs(ad_abeta(cerebrum_rois)));

        figstr = [figuresavename '-emciab'];
    yy =1*ballradius*emci_abeta;
    yy(Striatal_inds) = min(yy);
    [h, tiffname2] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-lmciab'];
    yy =1*ballradius*lmci_abeta;
    yy(Striatal_inds) = min(yy);
    [h, tiffname4] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-adab'];
    yy =1*ballradius*ad_abeta;
    yy(Striatal_inds) = min(yy);
    [h, tiffname6] = glassbrain_display(yy(cerebrum_rois), figstr);

figure; 
cnt = 0;
for i = 1:length(tsamples)
    subplot(length(tsamples),2,cnt+1); subimage(imread(tiffnamef{i}));  title(['A\beta spread, t = ' num2str(tsamples(i))]); axis off
    cnt = cnt+2;
end

subplot(length(tsamples),2,2); plot(tsol, Rall_emci, 'g-', tsol, Rall_lmci, 'r-', tsol, Rall_ad, 'b-');
title('Similarity index between model and ADNI patients AV45-PET (green: EMCI, red: LMCI, blue: AD)'); xlabel('Model time (a.u.)'); ylabel('Similarity index');
subplot(length(tsamples),2,4); subimage(imread(tiffname2));  title(['ADNI EMCI AV45 uptake, Rmax = ' num2str(R_emci)]); axis off
subplot(length(tsamples),2,6); subimage(imread(tiffname4));  title(['ADNI LMCI AV45 uptake, Rmax = ' num2str(R_lmci)]); axis off
subplot(length(tsamples),2,8); subimage(imread(tiffname6));  title(['ADNI AD AV45 uptake, Rmax = ' num2str(R_ad)]); axis off
 
% Second, ODE dynamics of xsol and empirical atrophy, tau
for i = 1:length(tsamples)
    figstr = [figuresavename '-ODEx' num2str(i)];
    [dummy, ind] = min(abs(tsol-tsamples(i)));
    yy = 3*ballradius*xsol(:,ind); % 0.1
%     yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamex{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
end
[Rall_emci, R_emci] = mycorr(xsol(cerebrum_rois,:), emci_atrophy(cerebrum_rois));
% [~, ind_emci] = max(Rall_emci);
% [R_emci] = mycorr(xsol(cerebrum_rois, ind_emci), abs(emci_atrophy(cerebrum_rois)));
[Rall_lmci, R_lmci] = mycorr(xsol(cerebrum_rois,:), lmci_atrophy(cerebrum_rois));
% [~, ind_lmci] = max(Rall_lmci);
% [R_lmci] = mycorr(xsol(cerebrum_rois, ind_lmci), abs(lmci_atrophy(cerebrum_rois)));
[Rall_ad, R_ad] = mycorr(xsol(cerebrum_rois,:), ad_atrophy(cerebrum_rois));
% [~, ind_ad] = max(Rall_ad);
% [R_ad] = mycorr(xsol(cerebrum_rois, ind_ad), abs(ad_atrophy(cerebrum_rois)));

[Rall_korea_namci, R_korea_namci] = mycorr(xsol(ii,:), emci_tau(ii));
% [~, ind_korea_namci] = max(Rall_korea_namci);
% [R_korea_namci] = mycorr(xsol(ii, ind_korea_namci), abs(emci_tau(ii)));
[Rall_korea_amci, R_korea_amci] = mycorr(xsol(ii,:), lmci_tau(ii));
% [~, ind_korea_amci] = max(Rall_korea_amci);
% [R_korea_amci] = mycorr(xsol(ii, ind_korea_amci), abs(lmci_tau(ii)));
[Rall_korea_ad, R_korea_ad] = mycorr(xsol(ii,:), ad_tau(ii));
% [~, ind_korea_ad] = max(Rall_korea_ad);
% [R_korea_ad] = mycorr(xsol(ii, ind_korea_ad), abs(ad_tau(ii)));

        figstr = [figuresavename '-emciatrophy'];
    yy =1*ballradius*emci_atrophy; %5
    [h, tiffname2] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_namcitau'];
    yy =ballradius*emci_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname3] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-lmciatrophy'];
    yy =1*ballradius*lmci_atrophy; %5
    [h, tiffname4] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_amcitau'];
    yy =ballradius*lmci_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname5] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-adatrophy'];
    yy =1*ballradius*ad_atrophy; %5
    [h, tiffname6] = glassbrain_display(yy(cerebrum_rois), figstr);
        figstr = [figuresavename '-korea_adtau'];
    yy =ballradius*ad_tau;
    yy(Striatal_inds) = min(yy);
    [h, tiffname7] = glassbrain_display(yy(cerebrum_rois), figstr);

figure; 
cnt = 0;
for i = 1:length(tsamples)
    subplot(length(tsamples),3,cnt+1); subimage(imread(tiffnamex{i}));  title(['Tau spread, t = ' num2str(tsamples(i))]); axis off
    cnt = cnt+3;
end
cnt = 0;
subplot(length(tsamples),3,cnt+2); h = plot(tsol, Rall_korea_namci, 'g-', tsol, Rall_korea_amci, 'r-', tsol, Rall_korea_ad, 'b-');  [h.LineWidth] = deal(2); 
title('Similarity index between model and ADNI AV1451 uptake'); xlabel('Model time (a.u.)'); ylabel('Similarity');
legend('EMCI', 'LMCI', 'AD');
subplot(length(tsamples),3,cnt+3); h = plot(tsol, Rall_emci, 'g-', tsol, Rall_lmci, 'r-', tsol, Rall_ad, 'b-');  [h.LineWidth] = deal(2); 
title('Similarity index between model and ADNI atrophy'); xlabel('Model time (a.u.)'); ylabel('Similarity');
legend('EMCI', 'LMCI', 'AD');

cnt = 3;
subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname2));  title(['ADNI EMCI atrophy, Rmax = ' num2str(R_emci)]); axis off
subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname3));  title(['EMCI tau uptake, Rmax = ' num2str(R_korea_namci)]); axis off
cnt = cnt + 3;
subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname4));  title(['ADNI LMCI atrophy, Rmax = ' num2str(R_lmci)]); axis off
subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname5));  title(['LMCI tau uptake, Rmax = ' num2str(R_korea_amci)]); axis off
cnt = cnt + 3;
subplot(length(tsamples),3,cnt+3); subimage(imread(tiffname6));  title(['ADNI AD atrophy, Rmax = ' num2str(R_ad)]); axis off
subplot(length(tsamples),3,cnt+2); subimage(imread(tiffname7));  title(['AD tau uptake, Rmax = ' num2str(R_korea_ad)]); axis off

 %% Plots of global pathology burden
 for i = 1:length(tsol)
    q = xsol(:,i);
    globaltau(i) = sum(q);
    q = fsol(:,i);
    globalab(i) = sum(q);
 end
 for i = 1:length(tsolnoint)
    q = xsolnoint(:,i);
    globaltaunoint(i) = sum(q);
 end
 figure; h = plot(tsol, globalab, 'r-', tsol, globaltau, 'b-', tsolnoint, globaltaunoint, 'g-');
 legend('A\beta', 'tau', 'tau (no int)');
 title('Global burden of modeled A\beta, tau and tau (no interaction)'); [h.LineWidth] = deal(2);


%% % Glassbrains of ODE dyanmics
figure;      
for i = 1:length(tsamples)
    figstr = [figuresavename '-ODEx' num2str(i)];
    [dummy, ind] = min(abs(tsol-tsamples(i)));
    yy = 3*ballradius*xsol(:,ind);
    %yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamex{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
    figstr = [figuresavename '-ODEf' num2str(i+length(tsamples))];
    yy = 4*ballradius*fsol(:,ind);
    %yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamef{i+length(tsamples)}] = glassbrain_display(yy(cerebrum_rois), figstr);
    
    figstr = [figuresavename '-ODExnoint' num2str(i)];
    [dummy, ind] = min(abs(tsolnoint-tsamples(i)));
    yy = 3*ballradius*xsolnoint(:,ind);
    yy = min(yy, mean(yy) + 1*std(yy));
    %yy = max(0, yy - 1*std(yy));
    [h, tiffnamexnoint{i}] = glassbrain_display(yy(cerebrum_rois), figstr);
%     figstr = [figuresavename '-ODEfnoint' num2str(i+length(tsamples))];
%     yy = 10*ballradius*fsolnoint(:,ind);
%     yy = min(yy, mean(yy) + 1*std(yy));
%     %yy = max(0, yy - 1*std(yy));
%     [h, tiffnamefnoint{i+length(tsamples)}] = glassbrain_display(yy(cerebrum_rois), figstr);
end
figure; 
cnt = 0;
for i = 1:length(tsamples)
    subplot(length(tsamples),3,cnt+1); subimage(imread(tiffnamef{i+length(tsamples)}));   axis off
    subplot(length(tsamples),3,cnt+2); subimage(imread(tiffnamex{i}));   axis off
    subplot(length(tsamples),3,cnt+3); subimage(imread(tiffnamexnoint{i}));  axis off
    cnt = cnt+3;
end
subplot(length(tsamples),3,1); title('A\beta spread');
subplot(length(tsamples),3,2); title('tau spread');
subplot(length(tsamples),3,3); title('tau (no interaction) spread');
cnt = 0;
for i = 1:length(tsamples)
    subplot(length(tsamples),3,cnt+1); ylabel(['t = ' num2str(tsamples(i))]);
    cnt = cnt+3;
end


%% 10) ADDED 11/11/20: AIC and fitlm between various models

[tsol, xfsol] = ode45(@myode, trange, xf0, odeopts); 
% [tsol, xfsol] = ode45(@myode2, trange, xf0, odeopts); 
xsol = xfsol(:, 1:nroi); xsol = xsol.';
fsol = xfsol(:, nroi+1:2*nroi); fsol = fsol.';

morder = 3;  %change acc to which model 3 = no-int, 4 = 1-way, 5 = 2-way

[Rall] = mycorr(xsol(ii,:), emci_atrophy(ii));
[~, ind] = max(Rall);
mdl_emci = fitlm(xsol(ii, ind), abs(emci_atrophy(ii))),
mdl_emci.ModelCriterion,
logL = mdl_emci.LogLikelihood;
myAIC_emci = -2*logL + 2*morder,

[Rall] = mycorr(xsol(ii,:), lmci_atrophy(ii));
[~, ind] = max(Rall);
mdl_lmci = fitlm(xsol(ii, ind), abs(lmci_atrophy(ii))),
mdl_lmci.ModelCriterion,
logL = mdl_lmci.LogLikelihood;
myAIC_lmci = -2*logL + 2*morder,

[Rall] = mycorr(xsol(ii,:), ad_atrophy(ii));
[~, ind] = max(Rall);
mdl_ad = fitlm(xsol(ii, ind), abs(ad_atrophy(ii))),
mdl_ad.ModelCriterion,
logL = mdl_ad.LogLikelihood;
myAIC_ad = -2*logL + 2*morder,

        

%% 8) Node scrambling

if node_scrambling
    Rqsave = zeros(ntrials,1);
    % atrophy
    q = ad_atrophy(cerebrum_rois);
    Remci = mycorr(xsol(cerebrum_rois,:), emci_atrophy(cerebrum_rois));
    Rlmci = mycorr(xsol(cerebrum_rois,:), lmci_atrophy(cerebrum_rois));
    Rad = mycorr(xsol(cerebrum_rois,:), ad_atrophy(cerebrum_rois));
    figure;
    subplot(1,2,1); h = plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-'); [h.LineWidth] = deal(2);
    title('R, Random permutations of regional atrophy'); hold on;
    for i = 1:ntrials
        q = q(randperm(length(cerebrum_rois)));
        Rq = mycorr(xsol(cerebrum_rois,:), q);
    %     Rq = overlap_measure(xsol(cerebrum_rois,:), q, pthr); 
        if mod(i,5)==0, plot(tsol, Rq, 'g-'); end
        Rqsave(i) = max(Rq(:));
    end
    legend( 'EMCI', 'LMCI', 'AD', 'permuted AD');
    subplot(1,2,2);
    hist(Rqsave, 20, 'FaceColor', 'g'); hold on; 
    line([max(Rad), max(Rad)], [0, ntrials/10], 'Color', 'k', 'LineWidth', 2); 
    line([max(Remci), max(Remci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2); 
    line([max(Rlmci), max(Rlmci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2); 
     title('Peak R of regional atrophy and randomized permutations');
     legend('permuted AD', 'AD', 'EMCI', 'LMCI');

    % Do same for Tau dist
    Remci = mycorr(xsol(ii,:), emci_tau(ii));
    Rlmci = mycorr(xsol(ii,:), lmci_tau(ii));
    Rad = mycorr(xsol(ii,:), ad_tau(ii));
    Rqsave = zeros(ntrials,1);
    q = ad_tau(ii);
    figure;
    subplot(1,2,1); h = plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-'); [h.LineWidth] = deal(2);
    title('R, Random permutations of regional Tau'); hold on;
    for i = 1:ntrials
        q = q(randperm(length(ii)));
        Rq = mycorr(xsol(ii,:), q); 
    %     Rq = overlap_measure(xsol(ii,:), q, pthr); 
        if mod(i,5)==0, plot(tsol, Rq, 'g-'); end
        Rqsave(i) = max(Rq(:));
    end
    legend( 'EMCI', 'LMCI', 'AD', 'permuted AD');
    subplot(1,2,2);
    hist(Rqsave, 20, 'FaceColor', 'g'); hold on; 
    line([max(Rad), max(Rad)], [0, ntrials/10], 'Color', 'k', 'LineWidth', 2); 
    line([max(Remci), max(Remci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2); 
    line([max(Rlmci), max(Rlmci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2); 
     title('Peak R of regional Tau and randomized permutations');
     legend('permuted AD', 'AD', 'EMCI', 'LMCI');

    % Do same for Abeta dist
    Remci = mycorr(fsol(cerebrum_rois,:), emci_abeta(cerebrum_rois));
    Rlmci = mycorr(fsol(cerebrum_rois,:), lmci_abeta(cerebrum_rois));
    Rad = mycorr(fsol(cerebrum_rois,:), ad_abeta(cerebrum_rois));
    Rqsave = zeros(ntrials,1);
    q = ad_abeta(cerebrum_rois);
    figure;
    subplot(1,2,1); h = plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-'); [h.LineWidth] = deal(2);
    title('R, Random permutations of regional A\beta'); hold on;
    for i = 1:ntrials
        q = q(randperm(length(cerebrum_rois)));
        Rq = mycorr(fsol(cerebrum_rois,:), q); 
    %     Rq = overlap_measure(fsol(cerebrum_rois,:), q, pthr); 
        if mod(i,5)==0, plot(tsol, Rq, 'g-'); end
        Rqsave(i) = max(Rq(:));
    end
    legend( 'EMCI', 'LMCI', 'AD', 'permuted AD');
    subplot(1,2,2);
    hist(Rqsave, 20, 'FaceColor', 'g'); hold on; 
    line([max(Rad), max(Rad)], [0, ntrials/10], 'Color', 'k', 'LineWidth', 2); 
    line([max(Remci), max(Remci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2); 
    line([max(Rlmci), max(Rlmci)], [0, ntrials/10], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2); 
     title('Peak R of regional A\beta and randomized permutations');
     legend('permuted AD', 'AD', 'EMCI', 'LMCI');
end

if edge_scrambling
    xsol_true = xsol;
    Remci = mycorr(xsol_true(ii,:), emci_tau(ii));
    Rlmci = mycorr(xsol_true(ii,:), lmci_tau(ii));
    Rad = mycorr(xsol_true(ii,:), ad_tau(ii));
    Rqsave = zeros(ntrials,3);
    x0 = zeros(nroi,1);
    x0(seeds) = 1;
    x0 = x0/norm(x0); 
    f0 = nc_met.*nc_app; 
    f0 = f0/norm(f0);
    xf0 = [zeros(size(nc_met)); zeros(size(nc_met))];  
    for i = 1:ntrials        
        Crand = C(:);
        randinds = randperm(length(Crand));
        Crand = Crand(randinds);
        Crand = reshape(Crand, nroi, nroi);
        Crand = tril(Crand,-1) + triu(Crand.',1);
        rowdegree = (sum(Crand, 2)).';
        coldegree = sum(Crand, 1);
        Lret = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Crand; % * diag(1./(sqrt(coldegree)+eps)) ;

        [tsol, xfsol] = ode45(@myode, trange, xf0, odeopts);
%         [tsol, xfsol] = ode45(@myode2, trange, xf0, odeopts);
        xsol = xfsol(:, 1:nroi); xsol = xsol.';
        fsol = xfsol(:, nroi+1:2*nroi); fsol = fsol.';

        [tmp, Rqsave(i,1)] = mycorr(xsol(ii,:), emci_tau(ii)); 
        [tmp, Rqsave(i,2)] = mycorr(xsol(ii,:), lmci_tau(ii)); 
        [tmp, Rqsave(i,3)] = mycorr(xsol(ii,:), ad_tau(ii)); 
    end

%         plot(tsol, Remci, 'r:',tsol, Rlmci, 'r-',  tsol, Rad, 'k-', tsol, R21, 'b:', tsol, R22, 'b-'); 
%     title('R, model versus empirical AD atrophy and its randomized node permutations'); hold on;
%     legend( 'EMCI', 'LMCI', 'AD', 'MCI- A\beta-', 'MCI- A\beta+', 'permuted AD');
    figure; 
    subplot(131); hist(Rqsave(:,1), 30, 'FaceColor', 'g'); xlabel('R_{max}', 'FontSize', 12); ylabel('Frequency', 'FontSize', 12); xlim([0.25 0.65]); ylim([0 550]); hold on; 
    line([max(Remci), max(Remci)], [0, ntrials/4], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2); 
     title('Peak R of randomized connectome permutations, EMCI', 'FontSize', 14);     hold off;
    subplot(132); hist(Rqsave(:,2), 30, 'FaceColor', 'g'); xlabel('R_{max}', 'FontSize', 12); ylabel('Frequency', 'FontSize', 12); xlim([0.25 0.65]); ylim([0 550]); hold on; 
    line([max(Rlmci), max(Rlmci)], [0, ntrials/4], 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2); 
     title('Peak R of randomized connectome permutations, LMCI');     hold off;
    subplot(133); hist(Rqsave(:,3), 30, 'FaceColor', 'g'); xlabel('R_{max}', 'FontSize', 12); ylabel('Frequency', 'FontSize', 12); xlim([0.25 0.65]); ylim([0 550]); hold on; 
    line([max(Rad), max(Rad)], [0, ntrials/4], 'Color', 'k', 'LineWidth', 2); 
     title('Peak R of randomized connectome permutations, AD');     hold off;
     
    %  restore correct C etc
    Cret = C; %retro
    rowdegree = (sum(Cret, 2)).';
    coldegree = sum(Cret, 1);
    Lret = eye(nroi) - diag(1./(sqrt(rowdegree.*coldegree)+eps)) * Cret; % * diag(1./(sqrt(coldegree)+eps)) ;

end


%% 9) ADDED NEW 11/20/19: repeat over all possible seed locations (single seeds)
if repeat_seeding
    
% Doing this for aMCI subjects only - the most clinically interesting group
% Keep this commneted out for routine runs - this takes FOREVER

    seedvec = [1:34, 69:77; 35:68, 78:86];
    f0 = nc_met.*nc_app; 
    f0 = f0/norm(f0);
    xf0 = [zeros(size(nc_met)); zeros(size(nc_met))];  
    odeopts = odeset('NonNegative', 1:length(xf0));

    for q1 = 1:size(seedvec, 2)      
        x0 = zeros(nroi,1);
        x0(seedvec(:,q1)) = 1;
        x0 = x0/norm(x0); 
        [tsol, xfsol] = ode45(@myode, trange, xf0, odeopts);
%         [tsol, xfsol] = ode45(@myode2, trange, xf0, odeopts);
        
        xsol = xfsol(:, 1:nroi); xsol = xsol.';
        fsol = xfsol(:, nroi+1:2*nroi); fsol = fsol.';
        [tmp, Rrepeatseeding(q1,1)] = mycorr(xsol(ii,:), emci_tau(ii)); 
        [tmp, Rrepeatseeding(q1,2)] = mycorr(xsol(ii,:), lmci_tau(ii)); 
        [tmp, Rrepeatseeding(q1,3)] = mycorr(xsol(ii,:), ad_tau(ii)); 
    end
save Rrepeatseeding.mat Rrepeatseeding;
figure; subplot(311); bar(Rrepeatseeding(:,1)); title('Repeat tau seeding: Rmax of EMCI');
subplot(312); bar(Rrepeatseeding(:,2)); title('Repeat tau seeding: Rmax of LMCI');
subplot(313); bar(Rrepeatseeding(:,3)); title('Repeat tau seeding: Rmax of AD');

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%  %%%%%%%%   INTERNAL FUNCTIONS  %%%%%%%%%%%

    function xfprime = myode(t, xf)
        lambda = 15; % was hand-tweaked to give desired time range of 20 years 
        % Remaining parameter values given below were obtained after grid
        % seach with 10 bins per parameter. To re-optimize on your data run
        % this function within nested for loops to create your own grid search
        beta = 1; %1 
        alpha = 2; %3
        gama = 0.2; %0.5
        
        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
        gfun = (t/lambda)*exp(-t/lambda);
%         gfun = 5/tmax;    % 5 for tau, 20 for amyloid
        xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0;
%         xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + eta*max(xu-xthr, 0); % + gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0 + eta*max(fu-fthr, 0) ;
        xfprime = xfprime(:);
    end

    function xfprime = myode2(t, xf)
        % This one has interactions going both ways
        lambda = 15; % was hand-tweaked to give desired time range of 20 years 
        % Remaining parameter values given below were obtained after grid
        % seach with 10 bins per parameter. To re-optimize on your data run
        % this function within nested for loops to create your own grid search

        beta = 1; %1 
        alpha = 3; %3
        gama = 0.1; %0.5
        
        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
        gfun = (t/lambda)*exp(-t/lambda);
%         gfun = 5/tmax;    % 5 for tau, 20 for amyloid
        xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0 + gama*(fu.*xu) ;
%         xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + eta*max(xu-xthr, 0); % + gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0 + eta*max(fu-fthr, 0) ;
        xfprime = xfprime(:);
    end


    function xfprime = myode_nointeraction(t, xf)
        % This one has no interactions between tau and amyloid
        lambda = 15; % 15
        beta = 1; %1
        alpha = 2; %*beta;  %3

        xu = xf(1:nroi);
        fu = xf(nroi+1:2*nroi);
        xu = xu(:); fu = fu(:);
        gfun = (t/lambda)*exp(-t/lambda);
        xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0; % + gama*(fu.*xu) ;
        xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0;
%         xfprime(1:nroi) =             -beta*Lret*xu  + alpha*gfun*x0 + eta*max(xu-xthr, 0) ; %+ gama*(fu.*xu) ;
%         xfprime(nroi+1:2*nroi) = -beta* Lret*fu + alpha*gfun*f0 + eta*max(fu-fthr, 0) ;
        xfprime = xfprime(:);
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