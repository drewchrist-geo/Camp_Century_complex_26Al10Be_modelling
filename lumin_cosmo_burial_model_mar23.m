clear
clearvars
clear figures

%% Luminescence refined 26Al/10Be burial modelling of Camp Century subglacial sediment
% This script models the burial and paleo-exposure history of the Camp
% Century subglacial sediment given the luminescence-constrained
% depositional age (406 +- 35 ka) of the upper most sediment (1059-4)
% during MIS 11.

% Using the observed 26Al and 10Be concentrations and uncertainties in the
% upper-most (1059-4) and the lower-most (1063-7) sediments, this script
% first performs a Monte Carlo simulation of nuclide concentrations and
% then calculates average nuclide concentrations weighted by measurement
% uncertainty.

% Then, the 26Al and 10Be concentrations of 1059-4 are corrected for 406 +-
% 35
% of burial. This yields the 26Al/10Be ratio of the upper sediment at
% the time of sediment deposition
% in a small surface stream, assuming the sediment was sufficiently buried
% (by either sediment or ice) to prevent additional nuclide production
% since 406 +- 35 ka.

% The next step calculates the inventory of nuclides that could accumulate
% during a range of paleo-exposure periods of the upper sediment given the
% 26Al/10Be ratio of the upper sediment at the time of deposition. This
% also constrains the maximum duration (~16 kyr) of ice-free surface
% exposure of the upper sediment when the pre-MIS11 inherited mean
% 26Al/10Be ratio is equal to 0.

% For the lower sediment (1063-7), this script corrects the 26Al/10Be ratio
% for 406 +- 35 kyr of burial PLUS an additional period of burial (0-16
% kyr) that is equivalent to the range of exposure durations of the upper
% sediment.

% For both sediments, the script propagates the uncertainties of the
% nuclide concentrations and the luminescence age.

% This script generates two plots: Figure 1 is a two-isotope plot (banana
% plot) of 10Be concentrations vs the 26Al/10Be ratios of the upper (blue)
% and lower sediment (red): 1) weighted mean observed values and 2)
% corrected for 406 +- 35 kyr of burial. For the upper sediment, the
% modelled 26Al/10Be and 10Be concs given the range of MIS 11
% paleo-exposure durations is shown (blues and greens). For the lower
% sediment, the 26Al/10Be and 10Be concs are shown given the additional
% period of burial (reds and yellows).

% Figure 2 contains two subplots: shows the paleo-exposure period of the
% uppper sediment vs 1) the inherited 26Al/10Be ratio, and 2)the inherited
% total burial history.




%% 1. Model variables
% This first section contains variables that can be adjusted to run the
% model such as the number of simulations, the IRSL depositional age, and
% the duration of paleo-exposure prior to the depositional age of the upper
% most sediment.
 
sim = 100; %number of simulations for weighted avg nuclide concentrations (observed) and the uncertainty of OSL burial ages. Change this as needed. 

irsl_bur = 406000; irsl_bur_err = 35000; % IRSL mean +/- 1 standard error, 
tbur_sim = irsl_bur_err.*randn(sim,1) + irsl_bur; %create vector of random IRSL ages using mean age and standard error.

% 1.1 Paleo-exposure scenarios
% Iterative calibration of t_exp showed that,
% given a surface production ratio of 7.3, paleo-exposure longer than
% 16,000 years results in negative mean inherited 26Al/10Be ratios prior to
% MIS 11 exposure. If the production ratio is 6.7, then paleo-exposure
% greater than 17,000 years result in negative mean inherited 26Al/10Be
% ratios prior to MIS 11 exposure.

t_exp = (0:1000:17000)' ;% Creates row of exposure periods from 0 to 17 kyr at 1 kyr increments. 
z_t = length(t_exp); %use this to set up blank matrices for burial calculations

% loading nuclide data and parameters
load cc_nuclides.mat %imports observed nuclide concentrations from sample 1059-4 & 1063-7  at Camp Century

% Nuclide production and decay parameters
P26 = 30.3; P10 = 4.15; %Al-26 and Be-10 production rates at SLHL (atoms/(g*yr)). 
% Greenland specific ratio 7.3: P26=30.3.  Global ratio 6.7: P26 = 27.8

tau26 = 1.01E6 ; tau10 = 2.02E6; %mean lives
lam26 = log(2)/717000; lam10 = log(2)/1387000; %decay rates

%% 2. UPPER SEDIMENT 1059-4 
% 2.1 Simulated observed nuclide concentrations In this section, we create
% n simulations of observed nuclide concentrations in sediment 1059-4 based
% on the concentration and uncertainty of each measurement (n=6). Then the
% mean nuclide concentration weighted by the inverse sqrt of each
% measurment's uncertainty.

% Monte Carlo simulation of nuclide concentrations in 1059-4, which
% incorporates the uncertainty of each measurement, creates a [obs sim]
% matrix
al_sim = al26_err_up.*randn(6,sim) + al26_up;
be_sim = be10_err_up.*randn(6,sim) + be10_up;
rat_sim = al_sim./be_sim; %26Al/10Be ratio. Not used in calculations but used for keeping track along the way.

% Weighting factors of upper sediment sample 1059-4
inverr_26 = ((sqrt(al26_err_up)).^-1); inverr_10 = ((sqrt(be10_err_up)).^-1); %winverse sqrt of measurement uncertainties
wf_26 = inverr_26./sum(inverr_26(:)); wf_10 = inverr_10./sum(inverr_10(:)); % vectors of weighting factors based on inverse sqrt

% Calculate weighted avg of each column that incorporates uncertainty of
% each simulated measurement (creates [1 sim] matrix)
al26_wmu = sum(al_sim.*wf_26,1)./sum(wf_26,1); al26_sd = std(al_sim,wf_26);%weighted average(wmu) and weighted sd of the 6 random values for each of the sims, needs to be propagated for error calcs
al26_obs_mu = mean(al26_wmu); al26_obs_sd = std(al26_wmu); %mean and sd of the 100 simulations: 3.3 +/- 0.07 x10^5 atoms

be10_wmu = sum(be_sim.*wf_10,1)./sum(wf_10,1); be10_sd = std(be_sim,wf_10); %weighted average (wmu) and weighted sd of the 6 random values for each of the sims, needs to be propagated for error calcs
be10_obs_mu = mean(be10_wmu); be10_obs_sd = std(be10_wmu); %mean and sd of the 100 simulations: 7.5 +/- 0.07 x10^4 atoms

rat_wmu = al26_wmu./be10_wmu; rat_sd = rat_wmu.*sqrt(((be10_sd./be10_wmu).^2)+((al26_sd./al26_wmu).^2)); %weighted avg (wmu) and sd of the 26Al/10Be ratio, based ratio uncertainty based on Al26 and Be10 simulated weighted avg and weighted sd
rat_obs_mu = mean(rat_wmu); rat_obs_sd = mean(rat_sd); % Total overall mean +/- sd nuclide ratio: ~4.4 +/- 0.6. Just for keeping track

t_bur = (-log(rat_wmu/7.3))/(lam26-lam10); %initial total burial history of each of the 100 sims
t_bur_mu = mean(t_bur); t_bur_sd = std(t_bur); % mean initial burial history: 1.1 +/- 0.1 x10^6 yrs

cov_obs = cov(be10_wmu,rat_wmu); %covariance of be10 and ratio weighted averages (2 x 2) matrix. Used to plot error ellipses

%% 2.2 Solve for nuclide concentrations and ratio at the time of sediment deposition at 390 ka assuming:
% 1. burial age = IRSL sediment depositional age (i.e. end of MIS 11)
% 2. upper sediment remained completely buried and no additional nuclide
% production has occured since deposition.


%Simulated nuclide concentrations and ratios based on IRSL age and
%uncertainty. This creates a [sim sim] matrix (the columns are the nuclide
%simulations from above, the rows are the prescribed burial scenarios [0:1:14] kyr
%applied to each set of nuclide variations

al26_0 = reshape((al26_wmu./exp(-tbur_sim/tau26)),[],1); al26_0sd = reshape((al26_sd./exp(-tbur_sim/tau26)),[],1); %correct the 26Al conc and uncertainty for irsl burial. flattens into a single vector [sim^2, 1] to simplify calcs later
be10_0 = reshape((be10_wmu./exp(-tbur_sim/tau10)),[],1); be10_0sd = reshape((be10_sd./exp(-tbur_sim/tau10)),[],1); %correct the conc for burial. flatten into a single vector [sim^2, 1] to simplify calcs later

al26_0mu = mean(al26_0); al26_0unc = mean(al26_0sd); %mean concs at time of sediment deposition
be10_0mu = mean(be10_0); be10_0unc = mean(be10_0sd); %mean concs at time of sediment deposition

rat_0 = al26_0./be10_0; rat_0unc = rat_0.*sqrt((be10_0sd./be10_0).^2+(al26_0sd./al26_0).^2); %nuclide ratio simulations that account for variability of both measurement variability and IRSL age uncertainty
rat_0mu = mean(rat_0(:)); rat_0sd = mean(rat_0unc(:)); %26Al/10Be = 5.4 +/- 0.8 at the time of sediment deposition at 406 +- 35 ka 

% total burial history at the time of sediment deposition at 406 +- 35 ka
t_bur_0 = (-log(rat_0/7.3))/(lam26-lam10); %total burial history of each simulation corrected for 406 +- 35 kyr burial
t_bur_0mu = mean(t_bur_0); t_bur_0sd = std(t_bur_0); % mean total burial history after corrected for 406 +- 35 kyr burial: 0.67 +/-0.06 Myr

cov_0 = cov(be10_0,rat_0); % covariance of irsl-corrected be10 and ratio weighted averages in (2 x 2) matrix, used for plotting error ellipse

al26_0m = repmat(al26_0,1,z_t); be10_0m = repmat(be10_0,1,z_t); % creates blank [sim*sim z_t] cell array of end MIS 11 concentrations to subtract nuclides from different exposure scenarios from

%% 2.3 Calculate inherted nuclide concs and ratios based on different paleoexposure scenarios during MIS 11
% Using the nuclide concs at the end of MIS 11 (before sediment deposition
% at 406 +- 35 ka), apply a variety of exposure scenarios (t_exp) during
% MIS 11. subscript: inh = inherited

% calculate a vector of nuclides that would accumulate under different
% exposure periods (t_exp)
al26exp = (P26/lam26)*(1-exp(-t_exp*lam26))'; be10exp = (P10/lam10)*(1-exp(-t_exp*lam10))'; 

% Loop to fill in [sim*sim z_t] cell array of different exposure scenarios
% to subtract from end MIS 11 nuclide concs
al26expm = ones(sim*sim,z_t).*al26exp; be10expm = ones(sim*sim,z_t).*be10exp;

% Create [sim sim z_t] cell arrays of initial inherited MIS 11
% concentrations and 26Al/10Be ratio (endMIS11 - paleoexposure MIS 11)
al26inh =  reshape((al26_0m - al26expm),[],z_t)'; %re arrange matrix such that all simulations for each exposure scenario are in one row.
al26inh_mu = mean(al26inh,2); al26inh_sd = std(al26inh,[],2); %mean +/- sd inherited 26Al concs before MIS 11 exposure

be10inh =  reshape(be10_0m - be10expm,[],z_t)'; %re arrange matrix such that all simulations for each exposure scenario are in one row.
be10inh_mu = mean(be10inh,2); be10inh_sd = std(be10inh,[],2); %mean +/- sd inherited 10Be concs before MIS 11 exposure

rat_inh = al26inh./be10inh; %inherited 26Al/10Be ratio
rat_inh_mu = mean(rat_inh,2); ratinh_unc = rat_inh_mu.*sqrt((be10inh_sd./be10inh_mu).^2+(al26inh_sd./al26inh_mu).^2); %vectors of mean +/- sd 26Al/10Be ratios of each exposure scenario
%for t_exp kyr the 26Al/10Be ratio decreases from 5.4 to 0

berat_inh_mu = horzcat(be10inh_mu,rat_inh_mu); %create matrix of be concs and 26Al/10Be ratios for plotting

covinhblk = ones(2,2,z_t); %covariance of all paleo-exposure scenarios. used for plotting error ellipses
for p = 1:size(covinhblk,3)
    covinh(:,:,p)=cov(be10inh(p,:),rat_inh(p,:));
end

% Calculate total burial history of sediment based on inherited nuclide
% concentrations (i.e. prior to MIS 11 exposure)

t_bur_inh = (-log(rat_inh/7.3))/(lam26-lam10); % array of burial ages of sediment before MIS 11 for all sims. This generates imaginary numbers for t_exp >13 kyr
t_bur_inh_mu = mean(t_bur_inh,2); t_bur_inh_sd = std(t_bur_inh,[],2);%vectors of mean +/- sd ratios of each exposure scenario

%% 3. LOWER SEDIMENT 1063-7
% In this section, we make similar claculations to those above for the
% lower sediment (weighted average observed nuclide concs, nuclide concs at
% timing of sediment deposition) but instead of paleo-exposure simulations
% we account for the additional decay of nuclides assuming the lower
% sediment remained buried or shielded from nuclide production during MIS
% 11. Throughout this section 'l' for 'lower' appended to variables.

% 3.1 Burial history corrected for up to 422 (406 + 16) kyr burial (assume this sediment was buried during MIS 11). 
tburl_sim = (tbur_sim)+t_exp'; %adds the exposure duration to the IRSL burial duration for a total burial duration to beginning of last exposure 

% Monte Carlo simulation of observed 26Al & 10Be concs from lower sediment based on individual measurement and uncertainty. Creates [3 sim] matrix
al_siml = al26_err_lo.*randn(3,sim) + al26_lo; be_siml = be10_err_lo.*randn(3,sim) + be10_lo; 
rat_siml = al_siml./be_siml;

%Weighting factors for lower sample 1063-7. 
inverr_26_l = (sqrt(al26_err_lo)).^-1; wf_26_l = inverr_26_l./sum(inverr_26_l(:));  %weighting factors (wf) by inverse sqrt of measurement uncertainty
inverr_10_l = (sqrt(be10_err_lo)).^-1; wf_10_l = inverr_10_l./sum(inverr_10_l(:));

% Weighted mean of each column that incorporates uncertainty of each measurement (creates [1 sim] matrix)
al26l_wmu = sum(al_siml.*wf_26_l)/sum(wf_26_l); al26l_sd = std(al_siml,wf_26_l);
al26l_obs_mu = mean(al26l_wmu); al26l_obs_sd = std(al26l_wmu);

be10l_wmu = sum(be_siml.*wf_10_l)/sum(wf_10_l); be10l_sd = std(be_siml,wf_10_l);
be10l_obs_mu = mean(be10l_wmu); be10l_obs_sd = std(be10l_wmu);

ratl_wmu = al26l_wmu./be10l_wmu; ratl_unc = ratl_wmu.*sqrt((be10l_sd./be10l_wmu).^2+(al26l_sd./al26l_wmu).^2);%weighted avg (wmu) of the 26Al/10Be ratio 
ratl_obs_mu = mean(ratl_wmu); ratl_obs_sd = mean(ratl_unc); % Total overall mean +/- sd nuclide ratio: ~1.8 +/- 0.4. Just for keeping track

t_burl = (-log(ratl_wmu/7.3))/(lam26-lam10); %total burial history of lower sediment
t_burl_mu = mean(t_burl); t_burl_sd = std(t_burl); %3.1 +/- 0.2 Myr

%stuff for plotting later
beratl_obs_mu = horzcat(be10l_obs_mu,ratl_obs_mu); %compiled mean observed 10Be and 26Al/10Be for plotting
covl_obs = cov(be10l_wmu,ratl_wmu); %used for plotting error ellipses

% 3.2 LOWER SEDIMENT inherited nuclide concentrations and burial history

%Simulated nuclide concentrations and ratios assuming that the lower
%sediment remained buried since 406+-35 ka and during the exposure of the
%upper sediment (t_exp).

%This creates a [sim sim] matrix (the columns are the nuclide variations
%from above, the rows are range of burial scenarios applied to each set of
%nuclide variations

al26l_0 = reshape((al26l_wmu./exp(-tbur_sim/tau26)),[],1); al26l_0sd = reshape((al26l_sd./exp(-tbur_sim/tau26)),[],1); %correct the 26Al conc and uncertainty for irsl burial. flattens into a single vector [sim^2, 1] to simplify calcs later
be10l_0 = reshape((be10l_wmu./exp(-tbur_sim/tau10)),[],1); be10l_0sd = reshape((be10l_sd./exp(-tbur_sim/tau10)),[],1); %correct the conc for burial. flatten into a single vector [sim^2, 1] to simplify calcs later

al26l_0mu = mean(al26l_0); al26l_0unc = mean(al26l_0sd); %mean conc and sd at time of depositon of upper sediment
be10l_0mu = mean(be10l_0); be10l_0unc = mean(be10l_0sd);

ratl_0 = al26l_0./be10l_0; ratl_0unc = rat_0.*sqrt((be10l_0sd./be10l_0).^2+(al26l_0sd./al26l_0).^2); %nuclide ratio simulations that account for variability of both measurement variability and IRSL age uncertainty
ratl_0mu = mean(ratl_0(:)); ratl_0sd = mean(ratl_0unc(:)); %2.2 +/- 1.2 at 406 ka when upper sediment was deposited.

%burial history of lower sediment at 406 +- 35 ka
t_burl_0 = (-log(ratl_0/7.3))/(lam26-lam10); %total burial history of each simulation corrected for 406 +- 35 kyr burial
t_burl_0mu = mean(t_burl_0); t_burl_0sd = std(t_burl_0); % mean total burial history after corrected for 406 +- 35 kyr burial: 2.7 +/-0.2 Myr

covl_0 = cov(be10l_0,ratl_0); %covariance of irsl-corrected be10 and ratio weighted averages in (2 x 2) matrix; used for plottong error ellipses

al26l_0m = repmat(al26l_0,1,z_t); be10l_0m = repmat(be10l_0,1,z_t); %creates blank [sim*sim z_t] cell array of end MIS 11 concentrations to subtract nuclides from different exposure scenarios from

% 3.3 Correcting for additional burial of lower sediment during MIS 11 
% calculate a vector of nuclides that would decay in the lower sediment under different burial
% periods while upper sed is exposed (t_exp)
al26l_bur = exp(-t_exp/tau26)'; be10l_bur = exp(-t_exp/tau10)'; 

% Loop to fill in [sim*sim z_t] cell array of different exposure scenarios to subtract from end MIS 11 nuclide concs
al26l_bur_m = ones(sim*sim,z_t).*al26l_bur; be10l_bur_m = ones(sim*sim,z_t).*be10l_bur;

%
% Create [sim sim z_t] cell arrays of initial inherited MIS 11
% concentrations and 26Al/10Be ratio (endMIS11 + nuclides that accumulate
% during additional burial)
al26l_inh =  reshape((al26l_0m./al26l_bur),[],z_t)'; %re arrange matrix such that all simulations for each exposure scenario are in one row.
al26l_inh_mu = mean(al26l_inh,2); al26l_inh_sd = std(al26l_inh,[],2); %mean +/- sd inherited 26Al concs before MIS 11 exposure

be10l_inh =  reshape(be10l_0m./be10l_bur,[],z_t)'; %re arrange matrix such that all simulations for each exposure scenario are in one row.
be10l_inh_mu = mean(be10l_inh,2); be10l_inh_sd = std(be10l_inh,[],2); %mean +/- sd inherited 10Be concs before MIS 11 exposure

ratl_inh = al26l_inh./be10l_inh; %inherited 26Al/10Be ratio
ratl_inh_mu = mean(ratl_inh,2); ratl_inh_unc = ratl_inh_mu.*sqrt((be10l_inh_sd./be10l_inh_mu).^2+(al26l_inh_sd./al26l_inh_mu).^2); %vectors of mean +/- sd 26Al/10Be ratios of each exposure scenario


% Calculate total burial history of sediment before MIS 11 exposure
% (inherited burial history)
t_burl_inh = (-log(ratl_inh/7.3))/(lam26-lam10); % array of burial ages of sediment before MIS 11 for all sims. This generates imaginary numbers for t_exp >13 kyr
t_burl_inh_mu = mean(t_burl_inh,2); t_burl_inh_sd = std(t_burl_inh,[],2);%vectors of mean +/- sd ratios of each exposure scenario


beratl_inh_mu = horzcat(be10l_inh_mu,ratl_inh_mu); %create matrix of be concs and 26Al/10Be ratios for plotting

covl_inhblk = ones(2,2,z_t); %covariance of all exposure scenarios, used for plotting error ellipses
for p = 1:size(covl_inhblk,3)
    covl_inh(:,:,p)=cov(be10l_inh(p,:),ratl_inh(p,:));
end
%
%% Compile outputs for plotting in next section

be_mu_all = vertcat(be10_obs_mu,be10_0mu,be10inh_mu,be10l_obs_mu,be10l_0mu,be10l_inh_mu); %all mean Be concs in obs, end MIS 11, inherited
be_err_all = vertcat(be10_obs_sd,be10_0unc,be10l_inh_sd,be10l_obs_sd,be10l_0unc,be10l_inh_sd);%all Be err in obs, end MIS 11, inherited
rat_mu_all = vertcat(rat_obs_mu,rat_0mu,rat_inh_mu,ratl_obs_mu,ratl_0mu,ratl_inh_mu); %all mean Al/Be ratios in obs, end MIS 11, inherited
rat_err_all = vertcat(rat_obs_sd,rat_0sd,ratinh_unc,ratl_obs_sd,ratl_0sd,ratl_inh_unc); %all mean Al/Be ratios in obs, end MIS 11, inherited

berat_mu_all=horzcat(be_mu_all,be_err_all,rat_mu_all,rat_err_all); %compiled Be and ratios for plotting


%% PLOTS

% Figure 1: banana plot showing the 10Be and 26Al/10Be ratios: 1) observed
% weighted means of 1059-4 (cool colors) and 1063-7 (warm colors), 2) timing of upper sediment
% deposition at 406 +- 35 ka IRSL corrected ratio and 10Be for both
% samples, 3) modelled inherited ratio and 10Be prior to exposure of upper
% sediment and additional burial of lower sediment (406 + 16 = 422 kyr).

%set up banana plot: burial isochrons
t = (0:1E3:4E6);
t_burial = (0:1:9).*1E6; %0.5 Myr increments
t_burma = (0:1:9)';
n26 = ((P26/lam26)*(1-exp(-lam26*t)))';
n10 = ((P10/lam10)*(1-exp(-lam10*t)))';

n26b = n26*exp(-(lam26*t_burial));
n10b = n10*exp(-(lam10*t_burial));
n2610b = (n26./n10).*exp(-((lam26-lam10)*t_burial));

bur_chron = (num2str(t_burma, '%.1f Myr burial'));

%color scheme
blue  = [44 62 80]/255; %1059-4 color 
red  = [231 76 60]/255 ;%1063-6 color
Cb = colormap(winter(length(t_exp)));
Cr = flip(colormap(autumn(length(t_exp))));

t_exp_kyr = t_exp/1000;
t_exp_leg =(num2str(t_exp, '%-d yr exp'));
clf

% plots of all the simulations
figure(1);

% 1059-4
hold on

scatter(be10_obs_mu,rat_obs_mu,25,'MarkerFaceColor',blue,'MarkerEdgeColor','k','DisplayName','Mean Measured 1059-4') %plots weighted mean observed ratio and 10Be concs.
hold on

for i=1:size(be10inh,1)%to plot mean simulated 10Be 26Al 1059-4 under different paleoexposure scenarios
        plot(be10inh_mu(i,:),rat_inh_mu(i,:),'o','MarkerFaceColor',Cb(i,:),'DisplayName',t_exp_leg(i,:),'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',1)
end
hold on

% plot 2-sigma error ellipses for observed, end MIS11, and inherited solutions
elpt_obs = ellipsedata(cov_obs,[berat_mu_all(1,1),berat_mu_all(1,3)],100,2); plot(elpt_obs(:,1:2:end),elpt_obs(:,2:2:end),'Color',[0.5 0.5 0.5]);
hold on

elpt_0 = ellipsedata(cov_0,berat_mu_all(:,2),100,2); plot(elpt_0(:,1:2:end),elpt_0(:,2:2:end),'Color',[0.5 0.5 0.5]);
hold on

 for r = 1:size(covinh,3) 
        elpt_inh(:,:,r) = ellipsedata(covinh(:,:,r),berat_inh_mu(r,:),100,2); 
 end
 for r = 1:size(elpt_inh,3)
         plot(elpt_inh(:,1:2:end,r),elpt_inh(:,2:2:end,r),'Color',[0.5 0.5 0.5]);
 end       
hold on

% 1063-7
scatter(be10l_obs_mu,ratl_obs_mu,25,'MarkerFaceColor',red,'MarkerEdgeColor','k','DisplayName','Measured 1063-7') %plots weighted mean observed ratio and 10Be concs.
hold on

%to plot mean simulated 10Be under different 390 +- 40 kyr burial + burial equiv to exposure duration in 1063-7
for k=1:size(be10l_inh,1) 
        plot(be10l_inh_mu(k,:),ratl_inh_mu(k,:),'-o','MarkerFaceColor',Cr(k,:),'DisplayName',t_exp_leg(k,:),'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',1)  
end
hold on 

elpt_obsl = ellipsedata(covl_obs,beratl_obs_mu,100,2); plot(elpt_obsl(:,1:2:end),elpt_obsl(:,2:2:end),'Color',[0.5 0.5 0.5])
hold on

 for s=1:size(covl_inh,3) 
      elpt_inhl(:,:,s) = ellipsedata(covl_inh(:,:,s),beratl_inh_mu(s,:),100,2);
 end
for s = 1:size(elpt_inhl,3)
    plot(elpt_inhl(:,1:2:end,s),elpt_inhl(:,2:2:end,s),'Color',[0.5 0.5 0.5])
 end
hold on

%banana plot isochrons
plot(n10,n2610b,'Color',[0.8 0.8 0.8]); %plots 1 myr burial isochrons
hold on
plot(n10,n2610b(:,1),'Color','k','DisplayName','Surface production'); %plots surface production ratio line at top

%formatting the plot
set(gca, 'YScale', 'log','XScale','log')
ylim([0.1 10]); xlim([1E4 1E5]);

xlabel('[10Be] (atoms/g)'); ylabel('26Al/10Be');
box on

set(gcf,'units','centimeters','position',[0.1,0.1,17.8,14])

%% Figure 2: MIS 11 paleo-exposure modelling

% plotting MIS 11 solution space of paleo-exposure duration vs. inherited
% 26Al/10Be ratio and inherited total burial history

x = t_exp; 
y = rat_inh_mu;  dy = ratinh_unc; 
y2 = t_bur_inh_mu'; dy2 = t_bur_inh_sd'; 
y3 = ratl_inh_mu'; dy3 = ratl_inh_unc'; 
y4 = t_burl_inh_mu;  dy4 = t_burl_inh_sd;

figure(2)
xlim_2 = [0 16000];
subplot(2,1,1); %inherited 26Al/10Be ratio vs paleo-exposure scenarios

p1 = plot(x,y,'-',x,y+(2*dy),'--',x,y-(2*dy),'--'); %sample 1059-4 plot of mean and uncertainty
set(p1,'Color',blue);
hold on

p2 = plot(x,y3,'-',x,y3+(2*dy3),'--',x,y3-(2*dy3),'--'); %sample 1063-7 plot of mean and uncertainty
set(p2,'Color',red);

ylabel('Inherited ^{26}Al/^{10}Be');
ylim([0 8])
xlim(xlim_2)
set(gca,'YTick',0:1:8,'Xtick',0:2000:18000);

subplot(2,1,2); % Plots inherited total burial history of the upper and lower sediment
p3 = plot(x,y2,'-',x,y2+(2*dy2),'--',x,y2-(2*dy2),'--'); %upper sediment 1059-4
set(p3,'Color',blue);
hold on

p4 = plot(x,y4,'-',x,y4+(2*dy4),'--',x,y4-(2*dy4),'--'); %lower sediment 1063-7
set(p4,'Color',red);
ylim([0 6E6])

xlabel('Surface exposure duration (yr) of 1059-4 prior to deposition');
set(gca,'YTick',0:1E6:6E6,'Xtick',0:2000:18000);
ylabel('Inherited total burial prior (Myr)');
xlim(xlim_2)
