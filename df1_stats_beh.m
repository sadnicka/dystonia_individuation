
%% All participants

%% |*DEMOGRAPHICS*|
clear all;
D = dload('subject_list_all.txt');
fprintf('N=%d musicians, N=%d dystonics\n',length(unique(D.name(D.group==1))), ... 
    length(unique(D.name(D.group==2))));
  
% |*AGE DIFFERENCES*|
x=pivottable(D.name,D.group,D.age,'nanmean');
fprintf('Mus: %2.3f,SD=%2.3f\n',nanmean(x(:,1)),nanstd(x(:,1)));
fprintf('Dys: %2.3f,SD=%2.3f\n\n',nanmean(x(:,2)),nanstd(x(:,2)));
ttest(x(:,1),x(:,2),2,'independent');    

% |*Tubiana-Champagne Score*|
x=pivottable(D.name,[],D.TCS,'nanmean','subset',D.group==2);
fprintf('Dys: %2.3f,SD=%2.3f\n',nanmean(x),nanstd(x));

% |*Time-since dystonia onset*|
x=pivottable(D.name,[],D.dystoniaOnset,'nanmean','subset',D.group==2);
fprintf('Dys: %2.3f,SD=%2.3f\n',nanmean(x),nanstd(x));


% Behavioural participants
% |*DEMOGRAPHICS*|
clear all;
D = load('ens_alldat_Peak.mat');
fprintf('N=%d musicians, N=%d dystonics\n',length(unique(D.subj(D.group==1))), ...                                           
    length(unique(D.subj(D.group==2))));

% |*TIME SINCE DYSTONIA ONSET*|
x=pivottable(D.subj,[],D.dystoniaOnset,'nanmean','subset',D.group==2);
fprintf('Mean=%2.4f, SD=%2.4f',mean(x),std(x));

% |*AGE DIFFERENCES*|
x=pivottable(D.subj,D.group,D.age,'nanmean');
fprintf('Mus: %2.3f,SD=%2.3f\n',nanmean(x(:,1)),nanstd(x(:,1)));
fprintf('Dys: %2.3f,SD=%2.3f\n\n',nanmean(x(:,2)),nanstd(x(:,2)));
ttest(x(:,1),x(:,2),2,'independent');                                       

%% Split half reliability (enslaving)
% |*Non-dystonic*|
clear all;
D = load('ens_splithalf_peakalldat.mat');
x = pivottable(D.SubjN,[],D.r,'mean(fisherz(x))','subset',D.group==1);
m = nanmean(x); 
SE= stderr(x); 
fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));

% |*dystonic*|
x = pivottable(D.SubjN,[],D.r,'mean(fisherz(x))','subset',D.group==2);
m = nanmean(x); 
SE= stderr(x); 
fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));


%% Split half reliability (mirroring)
% |*Non-dystonic*|
clear all;
D = load('mm_splithalf_peakalldat.mat');
x = pivottable(D.SubjN,[],D.r,'mean(fisherz(x))','subset',D.group==1);
m = nanmean(x); 
SE= stderr(x); 
fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));

% |*dystonic*|
x = pivottable(D.SubjN,[],D.r,'mean(fisherz(x))','subset',D.group==2);
m = nanmean(x); 
SE= stderr(x); 
fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));


%% Dystonia vs no dystonia (Enslaving)
% |*Left vs Left*|
clear all;
D = load('ens_alldat_Peak.mat');
x = pivottable(D.SubjN,D.group,mean(D.ens,2),'nanmean','subset',D.hand==1);
ttest(x(:,1),x(:,2),2,'independent');

% |*Right vs Right*|
x = pivottable(D.SubjN,D.group,mean(D.ens,2),'nanmean','subset',D.hand==2);
ttest(x(:,1),x(:,2),2,'independent');

% |*MIXED-ANOVA*|
clear all;
D = load('ens_alldat_Peak.mat');
anovaMixed(mean(D.ens,2),D.SubjN,'within',D.hand,{'hand'},...
                                'between',D.group,{'group'});

% |*STRENGTH DIFF*|
anovaMixed(mean(D.mvc,2),D.SubjN,'within',D.hand,{'hand'},...
                                'between',D.group,{'group'});
  

%% Dystonics vs non-dystonics (Mirroring)
% |*Left vs Left*|
clear all;
D = load('ens_alldat_Peak.mat');
x = pivottable(D.SubjN,D.group,mean(D.mmFull,2),'nanmean','subset',D.hand==1);
ttest(x(:,1),x(:,2),2,'independent');

% |*average mirroring*|
fprintf('Control: %2.4f\n',exp(nanmean(x(:,1))));
fprintf('Dystonia: %2.4f\n',exp(nanmean(x(:,2))));


% Dystonics vs non-dystonics (Mirroring)
% |*Right vs Right*|
x = pivottable(D.SubjN,D.group,mean(D.mmFull,2),'nanmean','subset',D.hand==2);
ttest(x(:,1),x(:,2),2,'independent');                            

% |*average mirroring*|
fprintf('Control: %2.4f\n',exp(nanmean(x(:,1))));
fprintf('Dystonia: %2.4f\n',exp(nanmean(x(:,2))));

%% Dystonia vs no dystonia (Enslaving pattern similarity)
% |*Left vs Left*|
clear all;
D = load('ens_patternsimilarity.mat');
x = pivottable(D.SubjN,D.group,D.r,'nanmean(fisherz(x))','subset',D.hand==1);
ttest(x(:,1),x(:,2),2,'paired');

% |*Right vs Right*|
x = pivottable(D.SubjN,D.group,D.r,'nanmean(fisherz(x))','subset',D.hand==2);
ttest(x(:,1),x(:,2),2,'paired');


%% No DYSTONIA %%
% |*ENSLAVED (LEFT VS RIGHT)*|
clear all;
D = load('ens_alldat_Peak.mat');
x=pivottable(D.SubjN,D.hand,mean(D.ens,2),'mean','subset',D.group==1);
ttest(x(:,1),x(:,2),2,'paired');

%% No DYSTONIA %%
% |*AVERAGE ENSLAVING*|
fprintf('Left: %2.4f\n',exp(mean(x(:,1))));
fprintf('Right: %2.4f\n',exp(mean(x(:,2))));

%% No DYSTONIA %%
% |*STRENGTH (LEFT VS RIGHT)*|
x=pivottable(D.SubjN,D.hand,mean(D.mvc,2),'mean','subset',D.group==1);
ttest(x(:,1),x(:,2),2,'paired');

%% No DYSTONIA %%
% |*AVERAGE STRENGTH*|
fprintf('Left: %2.4f\n',mean(x(:,1)));
fprintf('Right: %2.4f\n',mean(x(:,2)));


%% With DYSTONIA %%
% |*LEFT VS RIGHT*|
clear all;
D = load('ens_alldat_Peak.mat');
x=pivottable(D.SubjN,D.hand,mean(D.ens,2),'mean','subset',D.group==2);
ttest(x(:,1),x(:,2),2,'paired');

%% With DYSTONIA %%
% |*AVERAGE ENSLAVING*|
fprintf('Left: %2.4f\n',exp(mean(x(:,1))));
fprintf('Right: %2.4f\n',exp(mean(x(:,2))));

%% With DYSTONIA %%
% |*STRENGTH (LEFT VS RIGHT)*|
x=pivottable(D.SubjN,D.hand,mean(D.mvc,2),'mean','subset',D.group==2);
ttest(x(:,1),x(:,2),2,'paired');

%% With DYSTONIA %%
% |*AVERAGE STRENGTH*|
fprintf('Left: %2.4f\n',mean(x(:,1)));
fprintf('Right: %2.4f\n',mean(x(:,2)));

%% With DYSTONIA %%
% |*Correlation of enslaving with Tubiana-Champagne Score*|
x=pivottablerow(D.subj,[D.TCS, mean(D.ens,2)],'mean','subset',D.group==2);
lm = fitlm(x(:,1),x(:,2));

fprintf('Rsquared: %2.4f\n',sqrt(lm.Rsquared.Ordinary));
lm.disp;





