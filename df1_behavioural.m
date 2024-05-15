function varargout = df1_behavioural(what,varargin)
% Enslaving anlaysis for dystonia data set

rootDir         = [];
analysisDir     = [rootDir '/data_summaries'];


%% Experimental Parameters
ens_label       = {'1/2','1/3','1/4','1/5',...
                   '2/1','2/3','2/4','2/5',...
                   '3/1','3/2','3/4','3/5',...
                   '4/1','4/2','4/3','4/5',...
                   '5/1','5/2','5/3','5/4'};
mm_label        = {'1/1','1/2','1/3','1/4','1/5',...
                   '2/1','2/2','2/3','2/4','2/5',...
                   '3/1','3/2','3/3','3/4','3/5',...
                   '4/1','4/2','4/3','4/4','4/5',...
                   '5/1','5/2','5/3','5/4','5/5'};

handLabel       = {'Musician LH' , 'Musician RH', 'Dystonia LH', 'Dystonia RH'};

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
colours     = {[0,0,0],cb(4, :),[0 0 1],cb(1, :),[0 1 0]};
sty_grp     = colours([3,5]);
sty_control = colours(3);
sty_patient = colours([5]);



%% Analysis Functions
switch (what)
    
    case 'preprocessing'
    
        % run mov_correct for d07/d08; this corrected using HFB with LFB scale
        % location of function: /Volumes/ANNA/data/FingerPattern_dystonia/Individuation_EMG/data_raw
    case 'run_BEH' % basic behavioural routine, subset of orginal Individj2b_subj/_trial
        %%  BUG : unless directly in d08 data_raw folder recalibrated force will not be used
        D = dload(fullfile(indiviDir,'data','subject_list.txt')); % This has 18 subjects listed (not s05 or s07 'rapid' individual task  - not currently for inclusion)
        
        for i = 1:size(D.name,1);
            fprintf('%s\n',D.name{i});
            cd(fullfile(behDir,D.name{i}));
            Individjhu2b_subj_dystonia(D.name{i}, 0);
        end

        % subject 6, block 3, trial 7
    case 'excl_extraBN' % 2 sujects had extra block numbers these are deleted in this condition
    
        clear all
        D = load('IN2b_d05_ana.mat'); % poor execution of task BN=1
        d = getrow(D, D.BN~=1);
        d.BN = d.BN-1;
        save('IN2b_d05_ana_orig.mat', '-struct', 'D')
        save('IN2b_d05_ana.mat', '-struct', 'd')
        
        clear all
        D = load('IN2b_s08_ana.mat'); % Problems with EMG, BN=11 deleted
        d = getrow(D, D.BN~=11);
        save('IN2b_s08_ana_orig.mat', '-struct', 'D')
        save('IN2b_s08_ana.mat', '-struct', 'd') 
    case 'make_alldat'      %% make alldat structure for analysis
                            
        D           = dload(fullfile(indiviDir,'data','subject_list.txt'));
        D.dysFinger = [D.L1 D.L2 D.L3 D.L4 D.L5 D.R1 D.R2 D.R3 D.R4 D.R5];
        cd(analysisDir);
        
        S=[];
        for i = 1:size(D.name,1);
            s = D.name(i);
            Si=load(sprintf('IN2b_%s_ana.mat',s{1}));
            Si.SubjN        = repmat(i,size(Si.BN,1),1);
            Si.subj         = repmat(s,size(Si.BN,1),1);
            Si.group        = repmat(D.group(i), size(Si.BN,1), 1);            
            Si.instrum      = repmat(D.instrum(i), size(Si.BN,1), 1);
            Si.dysFinger    = repmat(D.dysFinger(i,:), size(Si.BN,1), 1);
            Si.isScanned    = repmat(D.scan(i), size(Si.BN,1), 1);
            Si.age          = repmat(D.age(i), size(Si.BN,1), 1);
            Si.gender       = repmat(D.gender(i), size(Si.BN,1), 1);
            Si.dystoniaOnset= repmat(D.dystoniaOnset(i), size(Si.BN,1), 1);
            Si.TCS          = repmat(D.TCS(i), size(Si.BN,1), 1);
                
                if size(Si.BN,1)>300, % check that 2 subjects with extra BN have 300trials
                    sprintf(s{1}) 
                end
                
             S = addstruct(S,Si);
        end;
        save(fullfile(analysisDir,'alldat.mat'),'-struct','S');
        varargout={S};
    case 'REL_splithalf'        %% Split-half reliability of enslaving/mirroring patterns 
        type        = 'peakENS';
        transform   = inline('log(x)');            
        splits      = {[1 3 5 7 9],[2 4 6 8 10]};
        vararginoptions(varargin,{'type','splits'});  

        if strcmp(type,'peakENS')
            ensType = 'ens';
            dim     = 'red';
        else
            ensType = 'mm';
            dim     = 'full';
        end;
        % 1. Load trial data 
        D           = load(fullfile(analysisDir,'alldat.mat'));
        
        % 2. Loop over every subject % hemisphere, estimate enslaving matrices for each
        % individual split
        S=[];
        for sn=unique(D.SubjN)'
            for h=1:2;
                Dh      = getrow(D,D.SubjN==sn & D.hand==h);

                D1      = df1_behavioural('ENS_estimateSlope',getrow(Dh,ismember(Dh.BN,splits{1})),'type',type);
                D2      = df1_behavioural('ENS_estimateSlope',getrow(Dh,ismember(Dh.BN,splits{2})),'type',type);

                split1   = transform(ens_squareform(D1.(ensType),'dim',dim));
                split2   = transform(ens_squareform(D2.(ensType),'dim',dim));
                r        = corr(split1',split2');

                Si.SubjN        = sn;
                Si.hand         = h;
                Si.group        = Dh.group(1);
                Si.split1       = split1;
                Si.split2       = split2;
                Si.r            = r;
                S               = addstruct(S,Si);    
            end;
        end;
        
        % 3. Make hand label for easy access and plotting
        S.handLabel       = zeros(length(S.SubjN),1);                                     
        idxC              = S.group == 1;
        S.handLabel(idxC) = S.hand(idxC);
        idxD              = S.group == 2;
        S.handLabel(idxD) = S.hand(idxD) + 2;
        
        varargout = {S};
        save(fullfile(analysisDir,sprintf('%s_splithalf_peakalldat.mat',ensType)),'-struct','S');
    case 'plot_splithalf'
        E = load('splithalf_ens.mat');
        M = load('splithalf_mm.mat');
        
        barplot(E.handLabel, E.r); ylabel('splithalf ens')
        ttest(E.r(E.handLabel==1), E.r(E.handLabel==2), 2, 'independent') % dom vs non dominant
        ttest(E.r(E.handLabel==3), E.r(E.handLabel==4), 2, 'independent') % non dystonic vs dystonic
    
        figure()
        barplot(M.handLabel, M.r); ylabel('splithalf mm')
        ttest(M.r(M.handLabel==1), M.r(M.handLabel==2), 2, 'independent') % dom vs non dominant
        ttest(M.r(M.handLabel==3), M.r(M.handLabel==4), 2, 'independent') % non dystonic vs dystonic        
        ttest(M.r(M.handLabel==2), M.r(M.handLabel==4), 2, 'independent') % non dystonic vs dystonic 
    case 'ENS_estimateSlope'                % Slope of enslaving vs applied force for each digit
        EValThresh  = 0.95;
        debug       = 0;
        D           = varargin{1};
        type        = 'peakENS';   % peakEns. enslaving, 
                                     % peakMM. mirror movements
        vararginoptions({varargin{2:end}},{'EValThresh','type','debug'});
        
        % 1. Turning on iteration limit warnings
        % Increasing iteration number in robustfit needs a system file edit in matlab
        if debug
            warning('on','stats:statrobustfit:IterationLimit')
        else
            warning('off','stats:statrobustfit:IterationLimit')
        end;
        
        % 2. Estimating strength of digit specific force/individuation
        % relationship over subjects/weeks
        T = [];
        
        for sn=unique(D.SubjN)'
            Ds = getrow(D,D.SubjN==sn);
            
            for h=unique(Ds.hand)'
                fprintf('S%d H%d\n',sn,h);                    
                Dh       = getrow(Ds, Ds.hand==h);
                
                % Determine whether enslaving or mirrored enslaving is
                % to be estimated
                
                if strcmp(type,'peakENS')
                    hand_idx = (h-1)*5+[1:5];
                else
                    hand_idx = (2-h)*5+[1:5];
                end;

                % Looping over each active digit
                
                for actFin=1:5
                    digit_idx       = Dh.digit == actFin;
                    f               = Dh.peakA(digit_idx);                            
                    meanDevP        = Dh.peakF(digit_idx,hand_idx);
                    meanDevPAct     = Dh.(type)(digit_idx);
                    f_sort          = sort(f);
                    clear ens ols rwls

                    % Active finger specific robust regression without intercept
                    try
                        [bAct]          = robustfit(f,meanDevPAct,[],[],'off');
                    catch err
                        disp(err);
                        bAct = nan;
                        break;
                    end
                    
                    % Passive finger specific robust regression without intercept
                    for passFin=1:5
                        [b,STATS]      = robustfit(f,meanDevP(:,passFin),[],[],'off');
                        ens(passFin)   = b;
                        ols(passFin)   = STATS.ols_s;
                        rwls(passFin)  = STATS.robust_s;

                        if (debug) 
                            subplot(2,3,passFin);
                            scatter(f,meanDevP(:,passFin),'ob','filled'); hold on;
                            plot(f_sort,f_sort*b,'r','LineWidth',1.5);  hold off;
                            title(['ActFin: ' num2str(actFin) ' PassFin: ' num2str(passFin)]);
                        end;
                    end;
                    if (debug)
                        subplot(2,3,6);
                        plot(1:5,ols,1:5,rwls);
                        legend('OLS','RWLS'); legend boxoff;
                        ylabel('RMSE');
                        % keyboard;
                    end;
                        
                    % Saving subject data
                    Ti.SubjN        = sn;
                    Ti.subj         = Dh.subj(1,:);
                    Ti.hand         = h;
                    Ti.digit        = actFin;
                    Ti.individ      = bAct;
                    Ti.group        = Dh.group(1);
                    Ti.isScanned    = Dh.isScanned(1);
                    Ti.dysFinger    = Dh.dysFinger(1,:);
                    Ti.instrum      = Dh.instrum(1);
                    Ti.dystoniaOnset= Dh.dystoniaOnset(1);
                    Ti.age          = Dh.age(1);
                    Ti.gender       = Dh.gender(1);
                    Ti.TCS          = Dh.TCS(1);
                    if strcmp(type,'peakENS')
                        Ti.ens  = ens;
                    else
                        Ti.mm   = ens;
                    end;
                    Ti.rlag             = pivottable([],Dh.targetForce,Dh.rlag,'mean');
                    
                    Ti.ols          = ols;
                    Ti.rwls         = rwls;
                    Ti.mvc          = Dh.mvc(1,(h-1)*5+actFin);

                    T = addstruct(T,Ti);                        
                end;
            end;
        end;
        varargout = {T};
        warning('on','stats:statrobustfit:IterationLimit');
    case 'ENS_make_alldat'                  %% Make an enslaving structure
        EValThresh  = 0.75;
        debug       = 0;
        vararginoptions({varargin{2:end}},{'EValThresh'});
        D = load(fullfile(analysisDir,'alldat.mat'));
%         D = getrow(D,D.rlag>=0.2);
        
        % Get enslaving for
        %   - same hand
        %   - opposite hand (mirror movements)
        D1 = df1_behavioural('ENS_estimateSlope',D,'type','peakENS',...
                                    'EValThresh',EValThresh,'debug',debug);
        D2 = df1_behavioural('ENS_estimateSlope',D,'type','peakMM',...
                                    'EValThresh',EValThresh,'debug',debug);        
        clear D;      
        D    = D1;
        D.mm = D2.mm;
        D.ensOverall    = D1.individ;
        D.mmOverall     = D2.individ;
        D               = rmfield(D,'individ');
        
        
        T = [];
        for sn=unique(D.SubjN)'
             Ds = getrow(D,D.SubjN==sn);
             for h=unique(Ds.hand)'
                Dh       = getrow(Ds,Ds.hand==h);
                if isempty(Dh)
                    break;
                end;
                if numel(Dh.ens)~=25
                    break
                end;
                
                % Saving subject data
                Ti.SubjN        = Dh.SubjN(1);
                Ti.subj         = Dh.subj(1);
                Ti.hand         = Dh.hand(1);
                Ti.group        = Dh.group(1);
                Ti.instrum      = Dh.instrum(1);
                Ti.isScanned    = Dh.isScanned(1);
                Ti.dysFinger    = Dh.dysFinger(1,:);
                Ti.dystoniaOnset= Dh.dystoniaOnset(1);
                Ti.age          = Dh.age(1);
                Ti.gender       = Dh.gender(1);
                Ti.TCS          = Dh.TCS(1);
                
                Ti.ensOverall   = mean(log(Dh.ensOverall));
                Ti.mmOverall    = mean(log(Dh.mmOverall));
                Ti.ensRaw       = ens_squareform(Dh.ens,'dim','red');
                Ti.ens          = log(ens_squareform(Dh.ens,'dim','red'));
                Ti.ensFull      = log(ens_squareform(Dh.ens,'dim','full'));
                Ti.mm           = log(ens_squareform(Dh.mm,'dim','red'));
                Ti.mmFull       = log(ens_squareform(Dh.mm,'dim','full'));
                Ti.mvc          = Dh.mvc';
                Ti.rlag         = mean(Dh.rlag);
                T = addstruct(T,Ti);                        
            end;
        end;        
        
        % Make hand label for easy plotting and access
        %   - 1. Control non dominant
        %   - 2. control dominant
        %   - 3. patient non paretic
        %   - 4. patient paretic 
        
        T.handLabel       = zeros(length(T.SubjN),1);                                     
        idxC              = T.group == 1;
        T.handLabel(idxC) = T.hand(idxC);
        idxP              = T.group == 2;
        T.handLabel(idxP) = T.hand(idxP) + 2;
        
        save(fullfile(analysisDir,'ens_alldat_Peak.mat'),'-struct','T');     
              
    case 'plot_meanmvc'
         D=load(fullfile(analysisDir,'ens_alldat.mat'));
         
         
%          figure;
%          subplot(1,3,[1 2]);
%          traceplot(1:5,D.mvc,'split',D.handLabel,'errorfcn','stderr','leg',handLabel);
%          set(gca,'XTick',1:5,'XLim',[0.8 5.2], 'fontSize', 12);
%          ylabel('MVC (N)'); xlabel('digit');
         
         % subplot(133);
         S = tapply(D,{'SubjN','group' 'hand'},{'mvc'});         
         figure;
         CAT.fillcolor={[1 1 1], [1 0 0.2], [1 1 1], [1 0 0.2]};
         myboxplot([S.group, S.hand], S.mvc, 'style_tukey', 'xtickoff') % linecolor', linecolor, 'fillcolor', fillcolor,'plotall',0, 'xtickoff');
         ylim([0 50]), set(gca,'YTick', [0:10:50]);ylabel('');
      
         % t-test to look for group differences
         S.avgMVC = mean(S.mvc,2);
         ttest(S.avgMVC(S.group==1),S.avgMVC(S.group==2),2,'independent');
    case 'STATS_OverallENS'
        % musician/dystonia dataset 
        cd(analysisDir)
        D=load(fullfile(analysisDir,'ens_alldat.mat'));
        
        % 1. Is there any difference in enslaving of active hand 
        % A) healthy musicians
        fprintf('ens musician LH vs RH\n')
        ttest   (   D.ensOverall(D.handLabel==1), ...
                    D.ensOverall(D.handLabel==2),  2, 'paired')
        % B) dystonia
        fprintf('ens dystonia LH vs RH\n')
        ttest   (	D.ensOverall(D.handLabel==3), ...
                    D.ensOverall(D.handLabel==4),  2, 'paired') 
                    
        % 2. Is there difference in mirroring 
         % A) healthy musicians
        fprintf('mm musician LH vs RH\n')
        ttest   (   D.mmOverall(D.handLabel==1), ...
                    D.mmOverall(D.handLabel==2),  2, 'paired')
        % B) dystonia
        fprintf('mm dystonia LH vs RH\n')
        ttest   (	D.mmOverall(D.handLabel==3), ...
                    D.mmOverall(D.handLabel==4),  2, 'paired') 
        
        % 3.  Take average across two hands for both groups
        D=tapply(D,{'SubjN','group', 'hand'},{'ens','mean(x,1)'},{'mmFull','mean(x,1)'});         
        
        figure;
             
        CAT.fillcolor = {[1 1 1],[1 0 0], [1 1 1], [1 0 0]};
        
        subplot(1,2,1), myboxplot([D.hand, D.group], mean(D.ens,2),'style_tukey', 'plotall',0, 'xtickoff'), ylim([-7 -2]), ylabel('')
        subplot(1,2,2), myboxplot([D.hand D.group], mean(D.mmFull,2), 'style_tukey', 'plotall',0, 'xtickoff'),  ylim([-7 -2]), ylabel('')
       
       

        
        
        S = [];
        for i=1:length(D.SubjN)
            d_i = getrow(D,i);
            Si.SubjN    = [i; i];
            Si.group    = [d_i.group; d_i.group];
            Si.ens      = [mean(d_i.ens); mean(d_i.mmFull)];
            Si.ensType  = [1; 2];
            S           = addstruct(S,Si);
        end;

        
        % ens compare groups 
        fprintf('ens\nmusician vs MD\n')
        ttest   (   S.ens   (S.group==1 & S.ensType==1), ...     
                    S.ens   (S.group==2 & S.ensType==1), 2, 'independent') 
        
        % mm compare groups
        fprintf('mm\nmusician vs MD\n')
        ttest   (   S.ens   (S.group==1 & S.ensType==2), ...     
                    S.ens   (S.group==2 & S.ensType==2), 2, 'independent') 
    case 'ENS_MMcorr'
        % 1. controls only (n =7)
        subplot(1,2,1)
        D = load(fullfile(analysisDir,'ens_alldat.mat'));
        T = getrow(D,D.group==1);        
        t=tapply(T,{'SubjN'},{'ensOverall','mean(x,1)'},{'mmOverall','mean(x,1)'}); 
        [r1, b1, t1, p1] = scatterplot(t.ensOverall, t.mmOverall, 'regression', 'linear', 'printcorr', 'markertype','o','markercolor',[0 0 1],'markersize',6,'intercept',0);
        
        % 2. dystonia only (n = 11)
        T = getrow(D,D.group==2);        
        t=tapply(T,{'SubjN'},{'ensOverall','mean(x,1)'},{'mmOverall','mean(x,1)'}); 
        hold on
        [r2, b2, t2, p2] = scatterplot(t.ensOverall, t.mmOverall, 'regression', 'linear', 'printcorr', 'markertype','o','markercolor',[1 0 0],'markersize',6,'intercept',0);
    
        %. all subjects
        subplot(1,2,2)
        d= tapply(D,{'SubjN'},{'ensOverall','mean(x,1)'},{'mmOverall','mean(x,1)'});
        [r3, b3, t3, p3] = scatterplot(d.ensOverall, d.mmOverall, 'regression', 'linear', 'printcorr','markersize',6,'intercept',0);

        keyboard
    case 'ENS_controlHandSimilarity'        % look at similarity of hands patterns in healthy participants
        dim = 10;
        vararginoptions(varargin,{'dim'});        
        
        % 1. Get data for non-paretic subjects only
        D = load(fullfile(analysisDir,'ens_alldat.mat'));
        T = getrow(D,D.group==1);
        
        % 2. Remove individual mean effects from either hand
        %       - svd for dim. reduction after mean removal
        T.ens   = bsxfun(@minus,T.ens,mean(T.ens,2));
%         [u,s,v] = svd(T.ens);
%         T.ens   = u(:,1:dim)*s(1:dim,1:dim);
%         MANOVA1(T.handLabel,T.ens);
        
        % 3. Basic enslaving pattern correlations
        %       - nondom with dom hand
        %       - nondom with nondom hand
        %       - dom with dom hand
        D = T;
        S = [];
       
        D1 = getrow(D,D.handLabel==1);
        D2 = getrow(D,D.handLabel==2);
            for m=1:2
                switch(m)
                    case 1
                        Si.r     = diag(corr(D1.ens',D2.ens'));
                    case 2
                        L = length(D1.SubjN);
                        r = zeros(L,1);
                        for l=1:L
                            r(l) = mean(corr(D1.ens(l,:)',D2.ens'));
                        end;
                        Si.r     = r;
                end;
                Si.SubjN = D1.SubjN;
                Si.model = m*ones(length(Si.r),1);
                S        = addstruct(S,Si);
            end;
       
        A           = S;
        S           = tapply(S,{'SubjN','model '},{'r','mean'});
        varargout   = {S,A};
        
        % Plot
        barplot([],S.r,'split',S.model);
        ylabel('Pearson r');
        set(gca,'XTickLabel',{'within','between'},...
                'YTick',0:0.1:0.8);
        set(gca,'FontSize',14);
            
        % Summary test
        %   - within bs between
        fprintf('\n\nSummary\n');
        fprintf('within vs between\n');
        ttest(fisherz(S.r(S.model==1)),fisherz(S.r(S.model==2)),2,'paired');
        keyboard
    case 'plot_ensShapeGroup'
        % musician/dystonia dataset 
        D=load(fullfile(analysisDir,'ens_alldat.mat'));
        D=tapply(D,{'SubjN','group', 'hand'},{'ens','mean(x,1)'},{'mmFull','mean(x,1)'});         

        figure;
        subplot(211);
        traceplot(1:20,D.ens,'split',[D.group, D.hand],'errorfcn','stderr','leg',{'musician','dystonic'});
        ylim([-6 -2])
        set(gca,'YTick',-8:-2, 'XTick',1:20,'XTickLabel',ens_label, 'FontSize', 12);
        ylabel('enslaving (a.u.)');
        
        subplot(212)
        traceplot(1:25,D.mmFull,'split',[D.group, D.hand],'errorfcn','stderr','leg',{'musician','dystonic'});
        ylim([-8 -4])
        set(gca,'YTick',-8:-2,'XTick',1:25,'XTickLabel',mm_label, 'FontSize', 12);
        ylabel('mirroring (a.u.)');
        
        dim = 6;
        [u,s,v]     = svd(D.ens);
        D.ensDim    = u(:,1:dim)*s(1:dim,1:dim);
        [u,s,v]     = svd(D.mmFull);
        D.mmDim     = u(:,1:dim)*s(1:dim,1:dim);
        
        D1 = getrow(D,ismember(D.group,[1 2])); % musician vs dystonics
        MANOVA1(D1.group,D1.ensDim);             % enslaving
        MANOVA1(D1.group,D1.mmDim);              % mirroring

%         D2 = getrow(D,ismember(D.group,[1 3])); % musician vs controls
%         MANOVA1(D2.group,D2.ensDim);             % enslaving
%         MANOVA1(D2.group,D2.mmDim);              % mirroring

        keyboard;
      
    case 'ENS_patternSimilarity'        % estimate similarity in patterns between controls and patients
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
       
        S = [];
        for hand=1:2
            Dc = getrow(D,D.group==1 & D.hand==hand);
            
            for i = unique(Dc.SubjN)'
                isC     = getrow(Dc,Dc.SubjN==i);
                notC    = getrow(Dc,Dc.SubjN~=i);

                % control with left out control
                Si.SubjN    = i;
                Si.r        = corr(isC.ens',mean(notC.ens,1)');
                Si.group    = 1;
                Si.hand     = hand;
                S           = addstruct(S,Si);

                % patients with left out control
                d = getrow(D,D.group==2 & D.hand==hand);
                Si.SubjN    = i;
                Si.r        = corr(isC.ens',mean(d.ens,1)');
                Si.group    = 2;
                Si.hand     = hand;
                S           = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(rootDir,'Individuation_EMG/analysis','ens_patternsimilarity.mat'),'-struct','S');
    case 'ENS_patternSimilarityIndividual'        % individual version of patient/control similarity analysis
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));

        S = [];
        for hand=1:2
            Dc = getrow(D,D.group==1 & D.hand==hand);
            
            for i = unique(Dc.SubjN)'
                isC     = getrow(Dc,Dc.SubjN==i);
                notC    = getrow(Dc,Dc.SubjN~=i);

                % control with left out control
                Si.SubjN    = i;
                Si.r        = corr(isC.ens',mean(notC.ens,1)');
                Si.group    = 1;
                Si.hand     = hand;
                S           = addstruct(S,Si);

                % patients with left out control
                d = getrow(D,D.group==2 & D.hand==hand);
                Si.SubjN    = d.SubjN;
                Si.r        = corr(d.ens',mean(notC.ens,1)');
                Si.group    = repmat(2,length(d.SubjN),1);
                Si.hand     = repmat(hand,length(d.SubjN),1);
                S           = addstruct(S,Si);
            end;
        end;
        T = tapply(S,{'SubjN','group','hand'},{'r','mean(fisherz(x))'});
        varargout = {T};
        % save(fullfile(rootDir,'Individuation_EMG/analysis','ens_patternsimilarity.mat'),'-struct','S');        
    case 'MM_patternSimilarityIndividual'        % individual version of patient/control similarity analysis
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));

        S = [];
        for hand=1:2
            Dc = getrow(D,D.group==1 & D.hand==hand);
            
            for i = unique(Dc.SubjN)'
                isC     = getrow(Dc,Dc.SubjN==i);
                notC    = getrow(Dc,Dc.SubjN~=i);

                % control with left out control
                Si.SubjN    = i;
                Si.r        = corr(isC.mmFull',mean(notC.mmFull,1)');
                Si.group    = 1;
                Si.hand     = hand;
                S           = addstruct(S,Si);

                % patients with left out control
                d = getrow(D,D.group==2 & D.hand==hand);
                Si.SubjN    = d.SubjN;
                Si.r        = corr(d.mmFull',mean(notC.mmFull,1)');
                Si.group    = repmat(2,length(d.SubjN),1);
                Si.hand     = repmat(hand,length(d.SubjN),1);
                S           = addstruct(S,Si);
            end;
        end;
        T = tapply(S,{'SubjN','group','hand'},{'r','mean(fisherz(x))'});
        varargout = {T};
        % save(fullfile(rootDir,'Individuation_EMG/analysis','ens_patternsimilarity.mat'),'-struct','S');        
    case 'ENS_MM_pattern_similarity'
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        D = getrow(D,ismember(D.handLabel,[2 4]));

        style.use('group_smallmarker');        

        % 0. plot trace plot for enslaving and mirroring patterns
        plt.subplot(2,3,[1 2]);
        plt.trace([],D.ens,'split',D.group);    % plot enslaving across groups
        plt.set(gca,'xtick',1:20,'ytick',[-5:1:-1],...
                    'yticklabel',round(exp([-5:1:-1]),3),'ratio','normal','xtick',1:20);
        plt.legend('northwest',{'healthy musicians (RH)','musicians'' with dystonia'});
        plt.labels([],'enslaved forces (N per 1N)','Enslaving pattern','A');
        ylim([-4.7 -2]);

        plt.subplot(2,3,[4 5]);
        plt.trace([],D.mmFull,'split',D.group);    % plot enslaving across groups
        plt.set(gca,'xtick',1:20,'ytick',[-7:1:-3],...
                    'yticklabel',round(exp([-7:1:-3]),3),'ratio','normal','xtick',1:25);
        plt.legend('northwest',{'healthy musicians (RH)','musicians'' with dystonia'});        
        plt.labels([],'mirrored forces (N per 1N)','Mirroring pattern','B');
        ylim([-6.6 -4]);

        % 1. plot correlation between patients and controls
        E = df1_behavioural('ENS_patternSimilarityIndividual');
        M = df1_behavioural('MM_patternSimilarityIndividual');
        E = getrow(E,E.hand==2);    % right hand moving
        M = getrow(M,M.hand==2);    % right hand moving

        E.type = ones(length(E.SubjN),1);
        M.type = ones(length(M.SubjN),1) + 1;
        S = addstruct(E,M);

        style.use('group_smallmarker');        

        plt.subplot(2,3,[3 6]);
        plt.line(S.type,S.r,'split',S.group,'plotfcn','fisherinv(mean(x))');
        ylim([0.3 0.9]);
        plt.set('ylim',[0.3 0.9],'ratio','normal','xticklabel',{'enslaving','mirroring'});
        plt.labels([],'Pearson''s r','Pattern similarity','C');        

        % 3. save figure
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');  

    case 'STATS_ens_pattern_similarity'
        % 1. get enslaving pattern correlations between patients and controls
        E = df1_behavioural('ENS_patternSimilarityIndividual');
        E = getrow(E,E.hand==2);    % right hand moving

        % average correlation between enslaving patterns for patients and controls
        x = pivottable(E.SubjN,E.group,E.r,'mean');
        fprintf('Enslaving pattern for patient correlation with control: %1.3f\n\n',fisherinv(nanmean(x(:,2))))

        % enslaving patterns for patients tested against controls
        fprintf('Patient versus control:\n');
        ttest(x(:,1),x(:,2),2,'independent');

    case 'STATS_mm_pattern_similarity'
        % 1. get mirroring pattern correlations between patients and controls
        M = df1_behavioural('MM_patternSimilarityIndividual');
        M = getrow(M,M.hand==2);    % right hand moving

        % average correlation between mirroring patterns for patients and controls
        x = pivottable(M.SubjN,M.group,M.r,'mean');

        % mirroring patterns for patients tested against controls
        fprintf('Patient versus control:\n');
        ttest(x(:,1),x(:,2),2,'independent');

    case 'ENS_dysfingers' % in dystonia looks for patterns across differnet combinations of dys / nondys fingers (symptomatic hand)
        
        D = load('ens_alldat_Peak.mat');
        
        d = getrow(D,D.hand==2 & D.group==2);    % right (2) symp hand, dystonia (2) group
        
        S = [];
        for i=1:length(d.SubjN)
            % per patients
            d_i  = getrow(d,i); % each individual with dystonia
            dys = find(d_i.dysFinger(6:10));
            non = find(~ismember(1:5,dys));
            x   = nonsym_squareform(d_i.ens,'dim','red');
            x   = x + diag(nan(5,1));
            y   = nonsym_squareform(d_i.mmFull,'dim','full');
            
           
            
            % combine in summary struct for plotting ens then mm
            Si.SN       = i;
           
            % diff combinatinos of dys and non fingers (ens mean)
            Si.fing_group      = 1;
            Si.ens      = nanmean(nanmean(x(dys,dys))); 
            S           = addstruct(S,Si);

            Si.fing_group = 2;
            Si.ens      = nanmean(nanmean(x(dys,non)));
            S           = addstruct(S,Si); 
            
            Si.fing_group  = 3;
            Si.ens      = nanmean(nanmean(x(non,dys)));
            S           = addstruct(S,Si);

            Si.fing_group = 4;
            Si.ens      = nanmean(nanmean(x(non,non)));
            S           = addstruct(S,Si);
            
            % dys and non fingers to mirroring (mean)
            Si.fing_group = 5;
            Si.ens      = nanmean(nanmean(y(dys, 1:5)));
            S           = addstruct(S,Si);

            Si.fing_group = 6;
            Si.ens      = nanmean(nanmean(y(non, 1:5)));
            S           = addstruct(S,Si);

        end;
        varargout = {S};
        
        % save(fullfile(dataDir,'ens_dysfingers.mat'),'-struct','S');

    case 'MM_dys' % look for evidence of mirror dystonia as per clinical definition
        
        D = load('ens_alldat_Peak.mat');
        
        % left hand (asympomatic) for ens / mm 
        d = getrow(D, D.hand==1 & D.group==2);

        S = [];
        for i=1:length(d.SubjN)
            % per patient
            d_i  = getrow(d,i); % index each individual with dystonia
            dys = find(d_i.dysFinger(6:10)); % index dystonic fingers
            non = find(~ismember(1:5,dys)); % index non
            x   = nonsym_squareform(d_i.ens,'dim','red');
            x   = x + diag(nan(5,1));
            y   = nonsym_squareform(d_i.mmFull,'dim','full');;
            
            % summary struct for plotting
            Si.SN       = i;
           
            Si.mm_group = 1;
            Si.mm       = nanmean(nanmean(y(1:5,dys)));
            S           = addstruct(S,Si);

            Si.mm_group = 2;
            Si.mm       = nanmean(nanmean(y(1:5,non)));
            S           = addstruct(S,Si); 
            
            
        end
            varargout = {S};
        
        % save(fullfile(dataDir,'mm_dysfingers.mat'),'-struct','S');
    
    case 'plot_ens_dysfingers'  
        
        clear all
        clf
        D = load('ens_alldat_Peak.mat'); % all data for group means
        r = load('ens_dysfingers.mat'); % finger spec analysis right instructed hand
        l = load('mm_dysfingers.mat'); % finger spec analysis left instructed hand
        

        %% FIGURE 2
        subplot(1,4,1) % enslaving
        right = getrow(D, D.hand==2);
        x_jitter = (rand(length(right.ensOverall),1)-0.5)*0.3;
        scatter(right.group+x_jitter, right.ensOverall, 'filled')
        hold on
        boxplot(right.ensOverall, right.group)
        ylim([-7 -1])

        subplot(1,4,2) %mirroring  
        scatter(right.group+x_jitter, right.mmOverall, 'filled')
        hold on
        boxplot(right.mmOverall, right.group)
        ylim([-7 -1])
        
        
        subplot(1,4,3:4)
        x_jitter = (rand(length(r.fing_group),1)-0.5)*0.3;
        scatter(r.fing_group+x_jitter, r.ens, 'filled')
        hold on
        boxplot(r.ens, r.fing_group)
        ylim([-7 -1])

        % STATS that accompany subplot 2:3
        %  pattern in symptomatic hand
        anova1([r.ens(r.fing_group==1), r.ens(r.fing_group==2), r.ens(r.fing_group==3), r.ens(r.fing_group==4)])
        spss = [r.ens(r.fing_group==1), r.ens(r.fing_group==2), r.ens(r.fing_group==3), r.ens(r.fing_group==4)];
        %  mirror pattern
        ttest(r.ens(r.fing_group==5), r.ens(r.fing_group==6), 2, 'paired')

        %% FIGURE 4
        figure()
        subplot(1,3,1) % enslaving
        left = getrow(D, D.hand==1);
        x_jitter = (rand(length(left.ensOverall),1)-0.5)*0.3;
        scatter(left.group+x_jitter, left.ensOverall, 'filled')
        hold on
        boxplot(left.ensOverall, left.group)
        ylim([-7 -1])
        

        subplot(1,3,2) %mirroring  
        scatter(left.group+x_jitter, left.mmOverall, 'filled')
        hold on
        boxplot(left.mmOverall, left.group)
        ylim([-7 -1])
      


        subplot(1,3,3)
        x_jitter = (rand(length(l.mm_group),1)-0.5)*0.3;
        scatter(l.mm_group+x_jitter, l.mm, 'filled')
        hold on
        boxplot(l.mm, l.mm_group)
        ylim([-7 -1])

        % STATS that accompany subplot 2:3
        %  pattern in symptomatic hand
        ttest(l.mm(l.mm_group==1), l.mm(l.mm_group==2), 2, 'paired')


    case 'get_matchSubj'     % get subjects who have measurements for both cortical and enslaving patterns
        region      = [1 11];   % bi-lateral S1
        stimtype    = 0;        % motor condition
        vararginoptions(varargin,{'region'});

        % 0. Select common participants which have both enslaving and
        % cortical distances        
        D = load(fullfile(rootDir,'RegionOfInterest','reg_distance_raw.mat'));
        D = getrow(D,ismember(D.region,region) & D.hand~=D.regSide & D.stimtype==stimtype); 
        S = load(fullfile(rootDir,'Individuation_EMG/analysis','ens_alldat.mat'));
        
        subj = intersect(S.subj,D.subj);
        
        T = [];
        for sn=1:length(subj)
            for h=1:2
                d_i = getrow(D,strcmp(D.subj,subj(sn)) & (D.hand+1)==h);
                Si = getrow(S,strcmp(S.subj,subj(sn)) & S.hand==h);
                
                ens = ens_squareform(Si.ens,'dim','red');
                ens = squareform((ens+ens')/2);
                mm  = ens_squareform(Si.mm,'dim','red');
                mm  = squareform((mm+mm')/2);

                Si      = rmfield(Si,{'ensRaw','ensFull','mmFull','mvc'});
                Si.reg  = d_i.region;
                Si.dist = d_i.dist;
                Si.ens  = ens;
                Si.mm   = mm;
                T       = addstruct(T,Si);
            end;
        end;
        save(fullfile(rootDir,'Individuation_EMG/analysis','matched_alldat.mat'),'-struct','T');
    case 'MODEL_corrBehImg'             % correlating behaviour (enslaving) and imaging
        D = load(fullfile(rootDir,'Individuation_EMG/analysis','matched_alldat.mat'));
        
        S       = [];
        rMean   = zeros(1,2);
        
        for g=[1:2]
            d_i = getrow(D,D.group==g);
            
            mEns = mean(d_i.ens,1);
            mMM = mean(d_i.mm,1);
            mDist = mean(d_i.dist,1);
            
            rMean(g)    = corr(mEns',mDist');
            
            d_i.ensMS    = bsxfun(@minus,d_i.ens,mEns);
            d_i.mmMS     = bsxfun(@minus,d_i.mm,mMM);
            d_i.distMS   = bsxfun(@minus,d_i.dist,mDist);
            
            d_i.r        = diag(corr(d_i.ensMS',d_i.distMS'));
            
            S = addstruct(S,d_i);
        end;
        
        CAT.markersize = 12;
        CAT.markertype = 'o';
        CAT.markercolor = {'b','r'};
        CAT.markerfill = {'b','r'};
        CAT.markercolor = {'b','r'};
        lineplot(S.group,S.r,'split',S.group,'markersize',20,'style_thickline','CAT',CAT);
        set(gca,'XLim',[0.5 2.5],'YLim',[-0.75 0.1]);
        ylabel('Pearsons r');
        set(gca,'XTickLabel',{'Musician','Dystonic'});
        
        
        line([0.5 1.5],[rMean(1) rMean(1)],'color','b','linewidth',1.5);
        line([1.5 2.5],[rMean(2) rMean(2)],'color','r','linewidth',1.5);
        keyboard;
    case 'Fig_clinicalTest'             % clinical measurement of finger involvement
        D           = dload(fullfile(indiviDir,'data','subject_list.txt'));
        instName    = {'piano','guitar'};
        
        % sort by group, then by instrument
        [~,idx] = sortrows([D.group D.instrum]);
        D       = getrow(D,idx);
        
        D.leftH = [D.L1 D.L2 D.L3 D.L4 D.L5];
        D.rightH= [D.R1 D.R2 D.R3 D.R4 D.R5];
        
        % make schedules for behavioural/fmri/tms        
        make_figure;
        h0  = subplot(1,6,1);
        imagesc_rectangle(D.group,'scale',[1 2],'YDir','reverse','MAP',colormap(gray));
        title('ID');
        h1  = subplot(1,6,2);
        imagesc_rectangle(D.instrum,'scale',[1 2],'YDir','reverse','MAP',colormap(hot));
        title('Instrument');
        h2  = subplot(1,6,3:4);
        imagesc_rectangle(D.leftH,'YDir','reverse','MAP',colormap(gray));
        title('left hand');
        h3  = subplot(1,6,5:6);
        imagesc_rectangle(D.rightH,'YDir','reverse','MAP',colormap(gray));
        title('right hand');
        
        set_graphics(gcf,'ytick',[],'xtick',[],'ax','normal','fontsize',12);
        set_graphics({h2,h3},'xtick',1:5,'xlabel','fingers','ax','normal','fontsize',10);
        set_graphics(h0,'ytick',1:length(D.name),'ax','normal','fontsize',10);
        set_graphics(h1,'ytick',1:length(D.name),'yticklabel',instName(D.instrum)','ax','normal','fontsize',10);

        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_2row');                      
              
    case 'FIG_forceTrace'

        D = load(fullfile(analysisDir,'ens_mm_force_trace.mat'));

        % cut data frame for when participant moved
        D.Force = D.Force(D.I2:D.I5,:);
        D.Force = bsxfun(@minus,D.Force,D.Force(1,:));        

        % force traces for left (active) and right (passive) hands
        %   dystonia d07 pressed with left thumb, at target force level of 0.75
        actF    = D.Force(:,1)';        
        passF   = D.Force(:,[5 9])';

        dt  = 0.05;
        t   = [1:1:length(actF)] * dt;

        % plot exemplary force traces for active, enslaving and mirrored forces
        
        style.use('forcetrace1');
        plt.subplot(141);
        plt.trace(t,actF,'leg',{'left thumb'},'leglocation','northwest');
        plt.labels('time (s)','Force (N)','Applied force','B');
        plt.set('ylim',[0 11]);


        style.use('forcetrace2');
        plt.subplot(142);
        plt.trace(t,passF,'split',[1;2],'leg',{'left little','right ring'},'leglocation','northwest');
        plt.labels('time (s)','Force (N)','Enslaved and mirrored forces','C');
        plt.set('ylim',[0 1.3]);

        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');        

    case 'FIG_logslope'
        D = load(fullfile(figureDir,'logslope.mat'));
        
        style.use('logslope');
        plt.subplot(131);
        plt.scatter(D.finst,D.fens,'intercept',0);
        plt.labels({'peak instructed','force (N)'},{'peak enslaved','force (N)'},[],'D');
        plt.set(gca,'xlim',[0 35],'ylim',[0 1.5],'xtick',0:10:30);        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'style','1x2');                
    case 'FIG_EnsDegree'
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        
        style.use('group_smallmarker');        
        plt.line(D.hand,mean(D.ens,2),'split',D.group,...
                                      'leglocation','southwest','leg',{'non-dystonic','dystonic'});

    case 'FIG_mvc'
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        
        plt.subplot(132);
        style.use('group');        
        plt.bar(D.group,mean(D.mvc,2),'split',D.hand,'leglocation','northwest');

        plt.set(gca,'xticklabel',{'LH','RH','LH','RH'},'ylim',[0 30]);
        plt.labels([],'strength (N)','MVF','B');
    case 'FIG_TCS'
        %D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        D = load('ens_alldat_Peak.mat')
        D = getrow(D,D.hand==2 & D.group==2);   % dystonia/right hand
      
        
        sty = style.custom('orange','markersize',6);
        plt.scatter(D.TCS,mean(D.ens,2),'style',sty);

        plt.set(gca,'ytick',-4:0.5:-2,...
                    'yticklabel',round(exp(-4:0.5:-2),3),'ax','square');
        plt.labels('TC score',{'enslaving','(N per 1N force)'},'TC score vs enslaving','C');      

    case 'FIG_patterns'                     % plot patterns for the different groups/hands
        %% 0. Get data
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        sc = [0 0.15];
        
        t = {'control (lh)','control (rh)','dys (lh)','dys (rh)'};
        h = plt.figure;
        for hl=1:4
            d_i = getrow(D,D.handLabel==hl);
            x = nonsym_squareform(mean(d_i.ens,1),'dim','red');
            
            plt.subplot(2,4,hl);
            plt.image(exp(x),'MAP', crameri('batlow'), 'scale',sc,'leg','none');
            title(t{hl});
        end;
        plt.set(h,'ax','square');
    case 'FIG_patternsLinePlot'            % plot patterns for the different groups/hands
        
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat')); % 36 rows with 2 rows per subject
        
        %1,2 = control L,R; 3,4-dystonia L,R;
        d = getrow(D,D.group==1);
        
        style.use('groupx3');
        plt.trace(1:20,D.ens,'split',D.handLabel,'leg',{'non-dystonic (RH)','dystonic (LH)','dystonic (RH)'},...
                                                 'leglocation','northwest');
        plt.set(gca,'xtick',1:1:20);
    case 'FIG_ensSimilarity'
        D = load(fullfile(analysisDir,'ens_patternsimilarity.mat'));
        
        style.use('group');        
        plt.line(D.hand,D.r,'split',D.group,'plotfcn','fisherinv(mean(fisherz(x)))',...
                                            'leglocation','south','leg',{'non-dystonic','dystonic'});

        plt.set(gca,'xticklabel',{'LH','RH'},'ylim',[0.4 0.9],'ytick',0.4:0.1:0.9,'ax','normal');
        plt.labels([],'Pearsons r',{'pattern','similarity'},'C');
    case 'FIG_MMDegree'
        D = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        
        style.use('group_smallmarker');        
        plt.line(D.hand,mean(D.mmFull,2),'split',D.group,...
                                      'leglocation','southwest','leg',{'non-dystonic','dystonic'});

        
    case 'PLOT_dysIndividRel'
        S = load(fullfile(analysisDir,'ens_dysindivid.mat'));
        
        style.use('altenmuller_replication');
        
        plt.figure;
        subplot(1,3,[1 2]);
        plt.box(S.grp,S.ens,'split',S.grp);
        
        plt.set(gca,'ylim',[-4.5 -2],'xtick',[]);
        plt.set(gca,'yticklabel',round(exp(get(gca,'ytick')),2));
        plt.labels([],'enslaving (N/1N)','enslaving split by movement type','D');
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'0.75x2');        
        
    case 'PLOT_EnslavingGroup'
        plt.subplot(121);
        df1_behavioural('FIG_EnsDegree');
        plt.set(gca,'xticklabel',{'LH','RH'},'ytick',[-5:0.5:-3],...
                                     'yticklabel',round(exp([-5:0.5:-3]),3),'ratio','normal','yprecision','%1.3f');
        
        plt.subplot(122);
        df1_behavioural('FIG_MMDegree');
        plt.set(gca,'xticklabel',{'LH','RH'},'ytick',[-5:0.5:-3],...
                                     'yticklabel',round(exp([-5:0.5:-3]),3),'ratio','normal','yprecision','%1.3f');

        plt.match('y');
        plt.labels([],{'forces in uninstructed fingers','(N per 1N force)'},'Enslaved','A',121);                                 
        plt.labels([],{'forces in uninstructed fingers','(N per 1N force)'},'Mirrored','B',122);                                 

        anot.hide_axis('y');
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');          
    case 'PLOT_EnslavingPatterns'
        df1_behavioural('FIG_patterns');
        plt.subplot(241);
        plt.panel('A');
        plt.subplot(2,4,5:7);
        df1_behavioural('FIG_patternsLinePlot');
        plt.panel('B');
        plt.subplot(2,4,8);
        df1_behavioural('FIG_ensSimilarity');
        plt.panel('C');
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'2x2');        
    case 'PLOT_Mirroring'
        df1_behavioural('FIG_MMDegree');
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x1');    

    case 'STATS_enslaving'
        % average enslaving across hands/group
        clear D
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.ens,2),'mean');
        x   = exp(nanmean(x,1));
        fprintf('Average enslaving across hands/group\n');
        fprintf('Control (LH): %1.3f\nControl (RH): %1.3f\nPatient (LH): %1.3f\nPatient (RH): %1.3f\n\n',x(1),x(2),x(3),x(4));

        % control vs patient (right hand)
        clear D x
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.ens,2),'mean');
        fprintf('control vs patient (RH)\n');
        ttest(x(:,2),x(:,4),2,'independent');
        fprintf('\n');

        % control vs patient (left hand)
        clear D x
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.ens,2),'mean');
        fprintf('control vs patient (LH)\n');
        ttest(x(:,1),x(:,3),2,'independent');
        fprintf('\n');

    case 'STATS_mirroring'
        % average mirroring across hands/group
        clear D
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.mmFull,2),'mean');
        x   = exp(nanmean(x,1));
        fprintf('Average mirroring across hands/group\n');
        fprintf('Control (LH): %1.3f\nControl (RH): %1.3f\nPatient (LH): %1.3f\nPatient (RH): %1.3f\n\n',x(1),x(2),x(3),x(4));

        % control vs patient (in left hand)
        clear D x
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.mmFull,2),'mean');
        fprintf('control vs patient (in LH)\n');
        ttest(x(:,1),x(:,3),2,'independent');
        fprintf('\n');

        % control vs patient (in right hand)
        clear D x
        D   = load(fullfile(analysisDir,'ens_alldat_Peak.mat'));
        x   = pivottable(D.SubjN,D.handLabel,mean(D.mmFull,2),'mean');
        fprintf('control vs patient (in RH)\n');
        ttest(x(:,2),x(:,4),2,'independent');
        fprintf('\n');        
    
    case 'STATS'        % compiles statistics for project to pdf
        
        % behavioural stats
        f   = which('df1_stats_beh.m');
        opt = struct('format','pdf','rootDir',d);
        publish(f,opt);
        
        % imaging stats
        f   = which('df1_stats_img.m');
        opt = struct('format','pdf','outputDir',d);
        publish(f,opt);
    otherwise
        fprintf('No such case...\n');
        
end

% Converts the enslaving structure from a matrix (5x5) into a vector and
% back to a matrix
function ens = ens_squareform(ens,varargin)
dim = 'full';   % full or red
vararginoptions(varargin,{'dim'});
nRow = size(ens);

switch(dim)
    case 'full'
        if nRow > 1
            ens = reshape(ens',1,numel(ens));
        else
            ens = reshape(ens,5,5)';
        end;
    case 'red'
        diag_idx = 1:6:25;
        if nRow > 1
            try
                ens = reshape(ens',1,numel(ens));
                ens(diag_idx) = [];
            catch
                ens(diag_idx) = [];
            end;
        else
            for i=diag_idx
                ens = [ens(1:i-1) 0 ens(i:end)];
            end;
            ens = reshape(ens,5,5)';
        end;
end
                    
% Estimates the cronbach's alpha
function alpha = cronbachsAlpha(x)
C       = cov(x);   
K       = size(x,2); 
varM    = trace(C)/K;
covM    = (C.*(1-eye(K))); 
covM    = sum(covM(:))/(K*(K-1));
alpha   = K*covM/(varM+(K-1)*covM);

% Estimates multiple t-tests on the specified columns of y
function ttestMultiple(y,tests,testLabels,varargin)
sided = 2;
type  = 'independent';
vararginoptions(varargin,{'sided','type'});
for i=1:length(tests)
    disp(testLabels{i});
    ttest(y(:,tests{i}(1)),y(:,tests{i}(2)),sided,type);
end;

