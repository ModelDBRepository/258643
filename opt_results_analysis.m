figs = 1;

load('output/030318_p3_WT_np140_ng8000_seed1.mat');

feat = cell2mat(feat);

all = [par err feat];

all(isnan(all))=0;
uni = unique(all,'rows');

%% Parameters for extracting information from the main matrix:
%       Which columns are parameters, errors, features?
%       What are the upper and lower parameter boundaries?
%       What is the 'acceptable' range for each feature?

num_params = 39;
num_errors = 11 ;
num_features = 10 ;

%par_lower = [ 1e-8 -60 1e-6 1e-5 1e-7 -30 0.5 1e-9 1e-9 -60 0.5 1e-6 1e-5 1e-7 -20 0.5 1e-8 1e-8 -25 0.5 1e-8 -100 -50 1e-8 1e-8 -50 0.5 1e-8 1e-8 -10 0.5 1e-8 1e-8 -45 0.6 1e-8 1e-8 -50 0.5 0.0001 0.0 0.001 1e-8 0.0001 ];
%par_upper = [ 1e-3 -50 1.0  10.0 0.1  10  1.5 1e-4 1e-4 -40 1.5 1.0  10.0 0.1  10  1.5 1e-3 1e-3 -5  1.5 1e-3 -70  -37 1e-3 1e-3 -30 1.5 0.01 0.01 10  1.5 1e-3 1e-3 -35 1.4 1e-3 1e-3 -30 1.5 0.001  2.0 0.1   1e-3 0.001  ];

% feats to use is a list of which features are important to know about, and
% the other ones will be supplementary features used for error calculation
feat_names = ['ISI90','vmvar','width','ir','vm','isi30','isi60','thresh','AP_amp','AHP_depth'];
feats_upper = [112.5986 4.0 3.8544 328.9248 -47.0248 349.3887 84.2387 -30.8887 73.8468 8.2918];
feats_lower = [51.5432  0.0 3.2613 253.9324 -53.9752 109.4347 56.0435 -35.2159 62.5831 4.2795];
%feat_names = ['oscfreq','oscamp','ohmicIR','sagratio'];
%feats_upper = [ 21.4 5.7 648 0.55 ];
%feats_lower = [ 13.4 1.7 296 0.45 ];

feats_range=abs(feats_upper-feats_lower);
feats_up_mid = feats_upper - (feats_range/4);
feats_low_mid = feats_lower + (feats_range/4);
feats_up_loose = feats_upper + (feats_range/4);
feats_low_loose = feats_lower - (feats_range/4);

%% Extract 'good' models
good = uni;
for i = 1:num_errors
    col = num_params+i;
    good = good( find (good(:,col) == 0 ) , : );
end

good = uni;
for i = 1:num_features
    col = num_params+num_errors+i;
    good = good( find( good(:,col) < feats_upper(i) ) , : );
    good = good( find( good(:,col) > feats_lower(i) ) , : );
end

very_good = uni;
for i = 1:num_features
    col = num_params+num_errors+i;
    very_good = very_good( find( very_good(:,col) < feats_up_mid(i) ) , : );
    very_good = very_good( find( very_good(:,col) > feats_low_mid(i) ) , : );
end

nearly_good = uni;
for i = 1:num_features
    col = num_params+num_errors+i;
    nearly_good = nearly_good( find( nearly_good(:,col) < feats_up_loose(i) ) , : );
    nearly_good = nearly_good( find( nearly_good(:,col) > feats_low_loose(i) ) , : );
end

bad = uni;

col = num_params+num_errors+2;
bad = bad( find( bad(:,col+1) > 0 | bad(:,col+2) > 0 ) , : );

%% Pars, errs, feats
uni_pars = uni(:,1:num_params);
uni_errs = uni(:,num_params+1:num_params+num_errors);
uni_feats = uni(:,num_params+num_errors+1:num_params+num_errors+num_features);

good_pars = good(:,1:num_params);
good_errs = good(:,num_params+1:num_params+num_errors);
good_feats = good(:,num_params+num_errors+1:num_params+num_errors+num_features);

bad_pars = bad(:,1:num_params);
bad_errs = bad(:,num_params+1:num_params+num_errors);
bad_feats = bad(:,num_params+num_errors+1:num_params+num_errors+num_features);

% Extract each generation
num_gens = 350;
pop_size = 144;
parents = {};
offspring = {};
%parents{1} = all(1:pop_size,:);
for i = 1:num_gens-1
    parents{i} = all(((((2*i)-1)*pop_size)+1):((((2*i)-1)*pop_size)+pop_size),:);
    offspring{i} = all((((2*i)*pop_size)+1):(((2*i)*pop_size)+pop_size),:);
end

%% Figure of feature space coverage
figure; subplot(3,3,1); scatter(good(:,27),good(:,28)); subplot(3,3,2); scatter(good(:,27),good(:,29)); subplot(3,3,3); scatter(good(:,27),good(:,30)); subplot(3,3,5); scatter(good(:,28),good(:,29)); subplot(3,3,6); scatter(good(:,28),good(:,30)); subplot(3,3,9); scatter(good(:,29),good(:,30));

%% PLSR
[Z_good, mu, sigma] = zscore(good_pars);
[Xloadings,Yloadings,Xscores,Yscores,betaPLS22,PLSPctVar] = plsregress(...
	((good_pars-mu)./sigma),good_feats(:,1:num_features),num_params);
yfitPLS = [ones(size(((good_pars-mu)./sigma),1),1) ((good_pars-mu)./sigma)]*betaPLS22;

% PLSR coefficients for each parameter, for each feature
figure; subplot(4,1,1); bar(betaPLS22(2:end,1)); subplot(4,1,2); bar(betaPLS22(2:end,2)); subplot(4,1,3); bar(betaPLS22(2:end,3)); subplot(4,1,4); bar(betaPLS22(2:end,4));

%% Figures of all population members in 2D feature space
if figs
    figure();
    subplot(1,2,1);
    scatter( uni(:,num_params+num_errors+feats_to_use(1)) , uni(:,num_params+num_errors+feats_to_use(2)) );
    title('all: amp vs freq');
    subplot(1,2,2);
    scatter( good(:,num_params+num_errors+feats_to_use(1)) , good(:,num_params+num_errors+feats_to_use(2)) );
    title('good: amp vs freq');
    saveas(gcf,'plots/feat_scatter.fig');
end

%% PCA

[wcoeff,score,latent,tsquared,explained]=pca(good_pars,'VariableWeights','variance');
c3=wcoeff(:,1:3);
coefforth = inv(diag(std(good_pars)))*wcoeff;
I=c3'*c3;
c3=coefforth(:,1:3);
I=c3'*c3;
cscores=zscore(good_pars)*coefforth;
par_names={'gbar_leak','e_leak','gbar_hcn','Vhalf_hcn','e_hcn','gbar_kv4_a','Vhalf_kv4_a','taumod_kv4_a','gbar_kerg','Vhalf_kerg','taumod_kerg','gbar_cat','Vhalf_cat','taumod_cat','gbar_cal','Vhalf_cal','taumod_cal','kf_cal','Pmax_cad','beta_cad','gbar_sk','km_sk'};
post_pca_rho = corr([score good_feats]);

%% Convert parents and offspring to PC score values

parents_score = {};
offspring_score = {};
%parents{1} = all(1:pop_size,:);
[Z_good,mu,sigma]=zscore(good_pars);
uni_score = bsxfun(@rdivide,bsxfun(@minus,uni(:,1:num_params),mu),sigma)*coefforth;
for i = 1:num_gens-1
    parents_score{i} = bsxfun(@rdivide,bsxfun(@minus,parents{i}(:,1:num_params),mu),sigma)*coefforth;
    offspring_score{i} = bsxfun(@rdivide,bsxfun(@minus,offspring{i}(:,1:num_params),mu),sigma)*coefforth;
end

parents_score_traj = zeros(num_gens-1,num_params);

% histograms of parameter values over time:
figure;
for i = 1:num_params
    subplot(num_params/2,2,i);
    mini = min(uni_pars(:,i));
    maxi = max(uni_pars(:,i));
    edges = mini:(maxi-mini)/100:maxi;
    for j=1:num_gens-1
        [N,edges] = histcounts(parents{j}(:,i),edges);
        parents_traj_hist{i}(:,j) = N;
    end
    imagesc(parents_traj_hist{i});
    set(gca,'YDir','normal');
end
figure;
for i = 1:num_params
    subplot(num_params/2,2,i);
    mini = min(uni_pars(:,i));
    maxi = max(uni_pars(:,i));
    edges = mini:(maxi-mini)/100:maxi;
    for j=1:num_gens-1
        [N,edges] = histcounts(offspring{j}(:,i),edges);
        offspring_traj_hist{i}(:,j) = N;
    end
    imagesc(offspring_traj_hist{i});
    set(gca,'YDir','normal');
end
% histograms of PC score values over time:
figure;
for i = 1:num_params
    subplot(num_params/2,2,i);
    mini = min(uni_score(:,i));
    maxi = max(uni_score(:,i));
    edges = mini:(maxi-mini)/100:maxi;
    for j=1:349
        [N,edges] = histcounts(parents_score{j}(:,i),edges);
        parents_score_traj_hist{i}(:,j) = N;
    end
    imagesc(parents_score_traj_hist{i});
    set(gca,'YDir','normal');
end
figure;
for i = 1:num_params
    subplot(num_params/2,2,i);
    mini = min(uni_score(:,i));
    maxi = max(uni_score(:,i));
    edges = mini:(maxi-mini)/100:maxi;
    for j=1:349
        [N,edges] = histcounts(offspring_score{j}(:,i),edges);
        offspring_score_traj_hist{i}(:,j) = N;
    end
    imagesc(offspring_score_traj_hist{i});
    set(gca,'YDir','normal');
end

%% Plot of progression of feature values over time
for i=1:num_gens-1
    temppar=parents{i};
    temppar(isnan(temppar))=0;
    amp_traj(i)=mean(temppar(:,26));
    amp_traj_min(i)=min(temppar(:,26));
    amp_traj_max(i)=max(temppar(:,26));
    freq_traj(i)=mean(temppar(:,27));
    freq_traj_min(i)=min(temppar(:,27));
    freq_traj_max(i)=max(temppar(:,27));
end

amp_traj(isnan(amp_traj))=0;
amp_traj_min(isnan(amp_traj_min))=0;
amp_traj_max(isnan(amp_traj_max))=0;
freq_traj(isnan(freq_traj))=0;
freq_traj_min(isnan(freq_traj_min))=0;
freq_traj_max(isnan(freq_traj_max))=0;

figure;
subplot(2,1,1);
h=fill([1 1:349 fliplr(1:349)],[amp_traj_min(end) amp_traj_max(1:349) fliplr(amp_traj_min(1:349))], 'r');
set(h,'facealpha',.5);
set(h,'EdgeColor','none');
hold on;
plot(amp_traj);
subplot(2,1,2);
h=fill([1 1:349 fliplr(1:349)],[freq_traj_min(end) freq_traj_max(1:349) fliplr(freq_traj_min(1:349))], 'r');
set(h,'facealpha',.5);
set(h,'EdgeColor','none');
hold on;
plot(freq_traj);

%% Scatterhist of All and Good results in PC space

groups(1:size(uni_score,1))={'All'};
groups(size(uni_score,1)+1:end)={'Good'};
figure;
scatterhist(sh_plot_dat(:,1),sh_plot_dat(:,2),'Group',groups,'Kernel','on');

%% Parallel histograms of parameters
num_bins = 50;
hist_mat_good_par = zeros(num_bins,num_params);
hist_mat_bad_par = zeros(num_bins,num_params);
for i = 1:num_params
    mini = min(uni_pars(:,i));
    maxi = max(uni_pars(:,i));
    edges = mini:(maxi-mini)/num_bins-1:maxi;
    [n,x] = hist(good_pars(:,i),edges);
    hist_mat_good_par(:,i) = n;
    [n,x] = hist(bad_pars(:,i),edges);
    hist_mat_bad_par(:,i) = n;
end
figure;
subplot(2,1,1);imagesc(hist_mat_good_par);set(gca,'YDir','normal');
subplot(2,1,2);imagesc(hist_mat_bad_par);set(gca,'YDir','normal');


%% Plot showing transformation from PC score into feature value:


%% Figures of PCA results
if figs
    figure();
    subplot(2,1,1)
    bar(post_pca_rho(num_params+feats_to_use(1),1:num_params));
    title('corr btwn PCs and amp');
    subplot(2,1,2)
    bar(post_pca_rho(num_params+feats_to_use(2),1:num_params));
    title('corr btwn PCs and freq');
    saveas(gcf,'plots/pc_corr_with_feats.fig');
end

% Which PCs are most correlated with the feats_to_use?
summed_post_pca_rho=abs(post_pca_rho(num_params+feats_to_use(1),1:num_params))+abs(post_pca_rho(num_params+feats_to_use(2),1:num_params));
summed_post_pca_rho_sorted=sort(summed_post_pca_rho);
top_3_pcs_inds=find(summed_post_pca_rho>=summed_post_pca_rho_sorted(end-2));
top_3_pcs_coefforth=coefforth(:,top_3_pcs_inds);
top_3_pcs_wcoeff=wcoeff(:,top_3_pcs_inds);


if figs
    figure();
    subplot(3,1,1)
    bar(top_3_pcs_coefforth(:,1));
    title(sprintf('Params in PC %d', top_3_pcs_inds(1)));
    subplot(3,1,2)
    bar(top_3_pcs_coefforth(:,2));
    title(sprintf('Params in PC %d', top_3_pcs_inds(2)));
    subplot(3,1,3)
    bar(top_3_pcs_coefforth(:,3));
    title(sprintf('Params in PC %d', top_3_pcs_inds(3)));
    saveas(gcf,'plots/params_in_top_3_pcs.fig');
end

% Figure of all 3 combinations of top 3 PCs correlated with features
% Feature 1:
if figs
    figure();
    for i=1:2
        for j=2:3
            if (i~=j)
                grid = gridize(score(:,top_3_pcs_inds(i)),score(:,top_3_pcs_inds(j)),good_feats(:,feats_to_use(1)),100);
                subplot(2,2,(j-3)+(2*i));
                imagesc(grid);
                set(gca,'YDir','normal');
                title(sprintf('PC %d vs PC %d, amplitude', top_3_pcs_inds(i), top_3_pcs_inds(j)));
                xlabel(sprintf('PC %d',top_3_pcs_inds(i)));
                ylabel(sprintf('PC %d',top_3_pcs_inds(j)));
            end
        end
    end
    saveas(gcf,'plots/amplitude_corr_with_top_PCs.fig');
end
    
% Feature 2:
if figs
    figure();
    for i=1:2
        for j=2:3
            if (i~=j)
                grid = gridize(score(:,top_3_pcs_inds(i)),score(:,top_3_pcs_inds(j)),good_feats(:,feats_to_use(2)),100);
                subplot(2,2,(j-3)+(2*i));
                imagesc(grid);
                set(gca,'YDir','normal');
                title(sprintf('PC %d vs PC %d, frequency', top_3_pcs_inds(i), top_3_pcs_inds(j)));
                xlabel(sprintf('PC %d',top_3_pcs_inds(i)));
                ylabel(sprintf('PC %d',top_3_pcs_inds(j)));
            end
        end
    end
    saveas(gcf,'plots/frequency_corr_with_top_PCs.fig');
end

%% Rotation of top PCs:

[rotated_top_3_pcs_coefforth,rotation_matrix_top_3_pcs] = rotatefactors(top_3_pcs_coefforth);
% cast PC values of each model from PC space into rotated component space:
X = zscore(good_pars);
rotated_scores = X*rotated_top_3_pcs_coefforth;

% Figure of parameters in 3 rotated components
if figs
    figure();
    subplot(3,1,1)
    bar(rotated_top_3_pcs_coefforth(:,1));
    title(sprintf('Params in rotated component 1'));
    subplot(3,1,2)
    bar(rotated_top_3_pcs_coefforth(:,2));
    title(sprintf('Params in rotated component 2'));
    subplot(3,1,3)
    bar(rotated_top_3_pcs_coefforth(:,3));
    title(sprintf('Params in rotated component 3'));
    saveas(gcf,'plots/params_in_rotated_components.fig');
end

% Figure of 3 rotated components vs features
% Feature 1:
if figs
    figure();
    for i=1:2
        for j=2:3
            if (i~=j)
                grid = gridize(rotated_scores(:,i),rotated_scores(:,j),good_feats(:,feats_to_use(1)),100);
                subplot(2,2,(j-3)+(2*i));
                imagesc(grid);
                set(gca,'YDir','normal');
                title(sprintf('Rotated component %d vs Rotated component %d, amplitude', i, j));
                xlabel(sprintf('Rotated component %d',i));
                ylabel(sprintf('Rotated component %d',j));
            end
        end
    end
    saveas(gcf,'plots/amplitude_corr_with_rotated_components.fig');
end

% Feature 2:
if figs
    figure();
    for i=1:2
        for j=2:3
            if (i~=j)
                grid = gridize(rotated_scores(:,i),rotated_scores(:,j),good_feats(:,feats_to_use(2)),100);
                subplot(2,2,(j-3)+(2*i));
                imagesc(grid);
                set(gca,'YDir','normal');
                title(sprintf('Rotated component %d vs Rotated component %d, frequency', i, j));
                xlabel(sprintf('Rotated component %d',i));
                ylabel(sprintf('Rotated component %d',j));
            end
        end
    end
    saveas(gcf,'plots/frequency_corr_with_rotated_components.fig');
end

%% using top 2 PCs only:
top_2_pcs_inds=find(summed_post_pca_rho>=summed_post_pca_rho_sorted(end-1));
top_2_pcs_coefforth=coefforth(:,top_2_pcs_inds);
top_2_pcs_wcoeff=wcoeff(:,top_2_pcs_inds);

[rotated_top_2_pcs_coefforth,rotation_matrix_top_2_pcs] = rotatefactors(top_2_pcs_coefforth);
% cast PC values of each model from PC space into rotated component space:
X2 = zscore(good_pars);
rotated_scores_2 = X2*rotated_top_2_pcs_coefforth;

% Figure of parameters in 3 rotated components
if figs
    figure();
    subplot(2,1,1)
    bar(rotated_top_2_pcs_coefforth(:,1));
    title(sprintf('Params in 2PC rotated component 1'));
    subplot(2,1,2)
    bar(rotated_top_2_pcs_coefforth(:,2));
    title(sprintf('Params in 2PC rotated component 2'));
    saveas(gcf,'plots/params_in_rotated_components.fig');
end

% Figure of 2 rotated components vs features
% Feature 1:
if figs
    figure();
    grid = gridize(rotated_scores_2(:,1),rotated_scores_2(:,2),good_feats(:,feats_to_use(1)),100);
    imagesc(grid);
    set(gca,'YDir','normal');
    title(sprintf('component %d vs component %d, amplitude', 1, 2));
    xlabel(sprintf('2PC Rotated component %d',1));
    ylabel(sprintf('2PC Rotated component %d',2));
    saveas(gcf,'plots/amplitude_corr_with_2PC_rotated_components.fig');
end

% Feature 2:
if figs
    figure();
    grid = gridize(rotated_scores_2(:,1),rotated_scores_2(:,2),good_feats(:,feats_to_use(2)),100);
    imagesc(grid);
    set(gca,'YDir','normal');
    title(sprintf('component %d vs component %d, frequency', 1, 2));
    xlabel(sprintf('2PC Rotated component %d',1));
    ylabel(sprintf('2PC Rotated component %d',2));
    saveas(gcf,'plots/frequency_corr_with_2PC_rotated_components.fig');
end


%% output data for plotting parameter range covered in the search

% for i = 1:16305
% for j = 1:22
% kerg_meta_params(i,j) = kerg_mean_par(j) + (good_kerg_mp(i,1)*wcoeff(j,20)) + (good_kerg_mp(i,2)*wcoeff(j,21)) + (good_kerg_mp(i,2)*wcoeff(j,22));
% end
% end

dlmwrite('plots/good_params.dat',good_pars,'\t');
%dlmwrite('kerg_meta_params.dat',kerg_meta_params,'\t');
min_good_params=min(good_pars);
max_good_params=max(good_pars);
range_good_params=abs(max_good_params-min_good_params);
range_pars = par_upper-par_lower;
good_params_norm = bsxfun(@rdivide,bsxfun(@minus,good_pars,par_lower),range_pars);
%kerg_meta_params_norm = bsxfun(@rdivide,bsxfun(@minus,kerg_meta_params,lb_par),range_par);
dlmwrite('plots/good_params_norm.dat',good_params_norm,'\t');
%dlmwrite('kerg_meta_params_norm.dat',kerg_meta_params_norm,'\t');

%% output parameter sets for running Neuron simulations of the models that 
%  are at the min, median and max values for the top 3 PCs

pc_params = [ ...
             good_pars(find(score(:,top_3_pcs_inds(1))==min(score(:,top_3_pcs_inds(1)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(1))==median(score(:,top_3_pcs_inds(1)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(1))==max(score(:,top_3_pcs_inds(1)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(2))==min(score(:,top_3_pcs_inds(2)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(2))==median(score(:,top_3_pcs_inds(2)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(2))==max(score(:,top_3_pcs_inds(2)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(3))==min(score(:,top_3_pcs_inds(3)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(3))==median(score(:,top_3_pcs_inds(3)))),:)' ...
             good_pars(find(score(:,top_3_pcs_inds(3))==max(score(:,top_3_pcs_inds(3)))),:)' ...
             ];
dlmwrite('pc_params.hoc',pc_params,',');

%% output parameter sets
out_params = [ ...
            good_pars(1,:)' ...
            good_pars(floor(size(good_pars,1)/8),:)' ...
            good_pars(floor((size(good_pars,1)/8)*2),:)' ...
            good_pars(floor((size(good_pars,1)/8)*3),:)' ...
            good_pars(floor((size(good_pars,1)/8)*4),:)' ...
            good_pars(floor((size(good_pars,1)/8)*5),:)' ...
            good_pars(floor((size(good_pars,1)/8)*6),:)' ...
            good_pars(floor((size(good_pars,1)/8)*7),:)' ...
            good_pars(floor((size(good_pars,1)/8)*8),:)' ...
            ];
dlmwrite('out_params.txt', out_params, ',');

%% save workspace
save('opt_results_workspace.mat');