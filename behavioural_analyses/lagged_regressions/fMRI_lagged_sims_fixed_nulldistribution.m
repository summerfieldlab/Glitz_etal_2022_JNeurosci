%NGarrett August 2020
%%simulations for gems (using model fits from RL model)

%close all, clear all, etc.
clear all; close all; clc;

%number of trials - 320 in fMRI, 480 behavioural
n_trials = 320;

%values market can take
trans_prob = [0.8 0.2];

%prob market change
prob_market_shuffle = 0.1;

%slope
beta_values = [2.323894972
1.128508163
1.548250133
2.084245655
0.732090315
2.496796767
2.326529244
2.389054206
1.295761545
2.104893838
1.171394482
2.219237407
0.803600627
2.639583119
1.504580593
0.628779231
2.116724153
2.223524316
0.754206548
2.264866912
1.88899658
2.571246807
1.99164536
1.989020614
0.972185873
2.905678456
1.482301597
1.444123078
2.185765793];

%learning rates (transformed)
eta_values = [0.876262561
0.297154443
0.20272044
0.234601369
0.373963012
0.841303508
0.718967862
0.812524927
0.248802866
0.794826209
0.788795485
0.486641756
0.793633779
0.886627821
0.644530069
0.552773005
0.659586278
0.529763453
0.680094387
0.56010234
0.661767874
0.683632161
0.777312167
0.515396825
0.534930672
0.873837491
0.229840169
0.38703844
0.87203604];

%wvals (transformed)
w_values = 0:1/28:1;

%number of subjects to run in each simulation
nsubs = length(w_values);

%number of sims
nsims = 500;

% %start measuring the time it takes to run the whole thing
% tic;

for sims = 1:nsims
    
    clearvars -except n_trials nsubs nsims trans_prob prob_market_shuffle eta_values beta_values w_values sims stats_dependent_fixed stats_independent_fixed stats_difference
    
    %print what sim you're on
    sims
    
    %loop over simultions
    for subs = 1:nsubs
        %% generating the data
        
        %read in the counterbalancing table from the actual experiment
        trialsCounterbalancing = readtable('TrialsCounterbalancingfMRI.txt');
        
        %pick parameters for this simulations from the set
        cur_sub = randperm(nsubs, 1);
        beta = beta_values(cur_sub);
        eta = eta_values(cur_sub);   
        w_raw = w_values(cur_sub);
        
        %select starting transition probabities; 1st entry is gems 1 and 2; 2nd
        %entry gems 3 and 4
        starting_prob = [trans_prob(randperm(length(trans_prob), 1)); trans_prob(randperm(length(trans_prob), 1)); trans_prob(randperm(length(trans_prob), 1))];
        
        %initalise probability estimates for each gem
        SR_gem = [0.5, 0.5, 0.5, 0.5];
        %initialise probability estimates for each
        SR_m = [0.5, 0.5];
        
        %read in and shuffle counterbalanced trials for the number of
        %blocks as done in experiment
        trialsCounterbalancingPart1 = trialsCounterbalancing(1:32,:);
        trialsCounterbalancingPart2 = trialsCounterbalancing(33:64,:);
        
        trialsCounterbalancingRun1a = trialsCounterbalancingPart1(randperm(size(trialsCounterbalancingPart1,1)),:);
        trialsCounterbalancingRun1b = trialsCounterbalancingPart2(randperm(size(trialsCounterbalancingPart2,1)),:);
        trialsCounterbalancingRun2a = trialsCounterbalancingPart1(randperm(size(trialsCounterbalancingPart1,1)),:);
        trialsCounterbalancingRun2b = trialsCounterbalancingPart2(randperm(size(trialsCounterbalancingPart2,1)),:);
        trialsCounterbalancingRun3a = trialsCounterbalancingPart1(randperm(size(trialsCounterbalancingPart1,1)),:);
        trialsCounterbalancingRun3b = trialsCounterbalancingPart2(randperm(size(trialsCounterbalancingPart2,1)),:);
        trialsCounterbalancingRun4a = trialsCounterbalancingPart1(randperm(size(trialsCounterbalancingPart1,1)),:);
        trialsCounterbalancingRun4b = trialsCounterbalancingPart2(randperm(size(trialsCounterbalancingPart2,1)),:);
        trialsCounterbalancingRun5a = trialsCounterbalancingPart1(randperm(size(trialsCounterbalancingPart1,1)),:);
        trialsCounterbalancingRun5b = trialsCounterbalancingPart2(randperm(size(trialsCounterbalancingPart2,1)),:);
        
        
        %set up the counterbalancing as used in the experiment - order of runs/blocks counterbalanced
        if rem(subs,2) == 1 %if uneven subject number
            trialsOverall = [trialsCounterbalancingRun1a; trialsCounterbalancingRun1b; trialsCounterbalancingRun2b;trialsCounterbalancingRun2a;...
                trialsCounterbalancingRun3a;trialsCounterbalancingRun3b;trialsCounterbalancingRun4b;trialsCounterbalancingRun4a; trialsCounterbalancingRun5a; trialsCounterbalancingRun5b]; %start Runs 1 and 3 with the dependent block
        else %if even subject numbers
            trialsOverall = [trialsCounterbalancingRun1b; trialsCounterbalancingRun1a; trialsCounterbalancingRun2a;trialsCounterbalancingRun2b;...
                trialsCounterbalancingRun3b;trialsCounterbalancingRun3a;trialsCounterbalancingRun4a;trialsCounterbalancingRun4b; trialsCounterbalancingRun5b; trialsCounterbalancingRun5a]; %start Runs 2 and 4 with the dependent block
        end
        
        %rename the read in variables
        forced_trial = trialsOverall.Var1;
        door_presented = trialsOverall.Var3;
        gem_presented = trialsOverall.Var2;
        correlated_block = trialsOverall.Var4;
        
        %populate which market is presented when (gems 1 and 2 = market 1; gems
        %3 and 4 = market 2)
        market_presented = gem_presented+10;
        market_presented(find(market_presented==11)) = 1;
        market_presented(find(market_presented==12)) = 1;
        market_presented(find(market_presented==13)) = 2;
        market_presented(find(market_presented==14)) = 3;
        
        Qmarket_presented = gem_presented+10;
        Qmarket_presented(find(Qmarket_presented==11)) = 1;
        Qmarket_presented(find(Qmarket_presented==12)) = 1;
        Qmarket_presented(find(Qmarket_presented==13)) = 2;
        Qmarket_presented(find(Qmarket_presented==14)) = 2;
        
        %outcome for s1 on each trial - same number of reward and loss per run
        rewards_short = repmat([1, -1], 1, height(trialsCounterbalancingRun1a)/2)';
        rewards = [];
        
        %shuffle this
        for i = 1:10
            rewards = [rewards; rewards_short(randperm(length(rewards_short)))];
        end
        
        %add outcome for s2 on each trial (always 0)
        rewards(:, 2) = repmat([0], 1, n_trials)';
        
        %initalise these as empty
        choose_left = []; EV_gem =[]; EV_otherGem = []; dv =[]; Qd = [];
        p_draw = []; prob_choose_left = []; state = []; outcome = [];
        
        w = w_raw;

        %loop over trials
        for trial = 1:n_trials
            
            %do not do this - just set w to w_raw above
            % adapt weights according to block
%             if correlated_block(trial) == 1 %dependent
%                 w = w_raw; 
%             else %independent
%                 w = 1-w_raw; 
%             end
            
            
            %change probabilites?
            if trial==1
                %three probabilities - 1 and 2 for gems in independent
                %block and 3 for dependent block
                probs(trial, 1:3) = starting_prob;
                
            else
                
                %set probabilites for this trial as the same as previous trial
                probs(trial, :) = probs(trial-1, :);
                
                %flip probability depending on random number generated
                if rand < prob_market_shuffle
                    probs(trial, market_presented(trial)) = trans_prob(randperm(length(trans_prob), 1));
                    while probs(trial-1, market_presented(trial))==probs(trial, market_presented(trial))
                        probs(trial, market_presented(trial)) = trans_prob(randperm(length(trans_prob), 1));
                    end
                    
                else
                end
                
            end
            
            %determine the other gem in a given context
            if gem_presented(trial) == 1
                other_gem_context = 2;
            elseif gem_presented(trial) == 2
                other_gem_context = 1;
            elseif gem_presented(trial) == 3
                other_gem_context = 4;
            elseif gem_presented(trial) == 4
                other_gem_context = 3;
            end
            
            if gem_presented(trial)<3
                
                index = 1;
                    
            else
                
                index = 2; 
                
            end
            
            %combined probability of getting to vault (from market and gem estimates)
            Vtot = w*SR_m(index) + (1-w)*SR_gem(gem_presented(trial));
            
            %expected value
            Qd(trial, 1:2) = [Vtot*rewards(trial, 1), (1-Vtot)*rewards(trial, 1)];
            
            %difference in value between doors
            dv(trial, 1) = Qd(trial, 1) - Qd(trial,2);
            
            %probability of choosing left (black)
            prob_choose_left(trial, 1) = 1 - (1/(1+exp(beta*dv(trial, 1))));
            
            %go ahead and implement l/r choice if the current trial is
            %not a forced trial
            if forced_trial(trial) == 0
                choose_left(trial, 1) = rand<=prob_choose_left(trial, 1);
            else
                if door_presented(trial) == 1 %black door
                    choose_left(trial,1) = 1;
                else
                    choose_left(trial,1) = 0; %white door
                end
            end
            
            %draw transition p
            p_draw(trial, 1) = rand<=probs(trial, market_presented(trial));
            
            %if choose left
            if choose_left(trial, 1)==1
                
                %get to s1 vs s2?
                if p_draw(trial, 1)
                    
                    %choose left and get to state 1
                    state(trial, 1) = 1;
                    
                    %ergo, increase Q estimate (estimate of p) for this gem
                    SR_gem(gem_presented(trial)) = SR_gem(gem_presented(trial)) + (1-w)*eta*(1-Vtot); %current gem/context
                    SR_m(index) = SR_m(index) + w*eta*(1-Vtot); %current context
                    
                    
                else
                    
                    %actually get to state 2 (despite choosing left)
                    state(trial, 1) = 2;
                    
                    %decrease Q estimate (estimate of p)
                    SR_gem(gem_presented(trial)) = SR_gem(gem_presented(trial)) + (1-w)*eta*(0-Vtot); %current gem/context
                    SR_m(index) = SR_m(index) + w*eta*(0-Vtot); %current context
                    

                end
                
                
            else %choose right
                
                %get to s1 vs s2?
                if p_draw(trial, 1)
                    
                    %choose right and get to state 2
                    state(trial, 1) = 2;
                    
                    %increase Q estimate (estimate of p)
                    SR_gem(gem_presented(trial)) = SR_gem(gem_presented(trial)) + (1-w)*eta*(1-Vtot); %current context
                    SR_m(index) = SR_m(index) + w*eta*(1-Vtot); %current context

                else
                    
                    %actually get to state 1 (despite choosing right)
                    state(trial, 1) = 1;
                    
                    %decrease Q estimate (estimate of p)
                    SR_gem(gem_presented(trial)) = SR_gem(gem_presented(trial)) + (1-w)*eta*(0-Vtot); %current context
                    SR_m(index) = SR_m(index) + w*eta*(0-Vtot); %current context
                   
                end
                
            end
            
            
            %trial outcome (1, -1 or 0 depending on state) - does not (should
            %not) influence update...
            outcome(trial, 1) = rewards(trial, state(trial, 1));
            
        if (SR_m(index)>1)
            SR_m(index) = 1;
        elseif (SR_m(index)<0)
            SR_m(index) = 0;
        end
        
        if (SR_gem(gem_presented(trial))>1)
            SR_gem(gem_presented(trial)) = 1;
        elseif (SR_gem(gem_presented(trial))<0)
            SR_gem(gem_presented(trial)) = 0;
        end
        
        end
        
        %% now we will actually analyse the generated data
       
        %set up empty variables
        MB_nback = ones(5,5)*NaN; %model-based information as to whether the dark door should be chosen
        same_gem_nback = ones(5,5)*NaN; %was the last trial's context identical to the current one
        diff_gem_nback = ones(5,5)*NaN; %or not

        %loop over trials and fill in these variables
        for i = 6:n_trials
            
            for n = 1:5
                
                same_gem_nback(i, n) = (gem_presented(i-n) == gem_presented(i));
                diff_gem_nback(i, n) = (gem_presented(i-n) ~= gem_presented(i));
                
                if choose_left(i-n,1) == 1 %if picked black on the last trial
                    if (state(i-n)==1) &  rewards(i,1) ==1 %state 1 outcome and reward trial (good because state 1 gain/loss)
                        MB_nback(i, n) = +1;
                    elseif (state(i-n)==2)  & rewards(i,1) == -1 %state 2 outcome and loss trial (good because state 2 zero outcome)
                        MB_nback(i, n) = +1;
                    elseif (state(i-n)==1) &rewards(i,1) == -1 %state 1 and loss trial (bad as state 1 is gain/loss)
                        MB_nback(i, n) = -1;
                    elseif (state(i-n)==2)  &rewards(i,1) == 1 %state 2 and gain trial (bad as state 2 is zero)
                        MB_nback(i, n) = -1;
                    end
                elseif choose_left(i-n,1) ==0 %picked white on last trial
                    if (state(i-n)==1) & rewards(i,1) ==1
                        MB_nback(i, n) = -1;
                    elseif (state(i-n)==2)  & rewards(i,1) == -1
                        MB_nback(i, n) = -1;
                    elseif (state(i-n)==1) & rewards(i,1) == -1
                        MB_nback(i, n) = +1;
                    elseif (state(i-n)==2) & rewards(i,1) == 1
                        MB_nback(i, n) = +1;
                    end
                else
                    MB_nback(i,n) = NaN;
                end
                
                %did you receive positive evidence for picking dark door
                positive_evidence_nback_1 = MB_nback(:, 1);
                positive_evidence_nback_2 = MB_nback(:, 2);
                positive_evidence_nback_3 = MB_nback(:, 3);
                positive_evidence_nback_4 = MB_nback(:, 4);
                positive_evidence_nback_5 = MB_nback(:, 5);
                
                %set up the same context variable
                MB_same_1_back = positive_evidence_nback_1;
                MB_same_2_back = positive_evidence_nback_2;
                MB_same_3_back = positive_evidence_nback_3;
                MB_same_4_back = positive_evidence_nback_4;
                MB_same_5_back = positive_evidence_nback_5;
                
                %replace those trials on which previous trial not same with
                %0
                MB_same_1_back(same_gem_nback(:, 1) == 0) =0;
                MB_same_2_back(same_gem_nback(:, 2) == 0) =0;
                MB_same_3_back(same_gem_nback(:, 3) == 0) =0;
                MB_same_4_back(same_gem_nback(:, 4) == 0)=0;
                MB_same_5_back(same_gem_nback(:, 5) == 0)=0;
                
                %make same-context into matrix
                MB_same_nback(i,1) = MB_same_1_back(i);
                MB_same_nback(i,2) = MB_same_2_back(i);
                MB_same_nback(i,3) = MB_same_3_back(i);
                MB_same_nback(i,4) = MB_same_4_back(i);
                MB_same_nback(i,5) = MB_same_5_back(i);
                
                %set up same variable for other context
                MB_diff_1_back = positive_evidence_nback_1;
                MB_diff_2_back = positive_evidence_nback_2;
                MB_diff_3_back = positive_evidence_nback_3;
                MB_diff_4_back = positive_evidence_nback_4;
                MB_diff_5_back = positive_evidence_nback_5;
                
                %replace those trials where previously same with 0
                MB_diff_1_back(same_gem_nback(:, 1) == 1) = 0;
                MB_diff_2_back(same_gem_nback(:, 2) == 1) = 0;
                MB_diff_3_back(same_gem_nback(:, 3) == 1) = 0;
                MB_diff_4_back(same_gem_nback(:, 4) == 1) = 0;
                MB_diff_5_back(same_gem_nback(:, 5) == 1) = 0 ;
               
                %make other-context into matrix
                MB_diff_nback(i,1) = MB_diff_1_back(i);
                MB_diff_nback(i,2) = MB_diff_2_back(i);
                MB_diff_nback(i,3) = MB_diff_3_back(i);
                MB_diff_nback(i,4) = MB_diff_4_back(i);
                MB_diff_nback(i,5) = MB_diff_5_back(i);
                
            end
            
        end
        
        %take the difference between the mean values of the same context's model
        %based value and different context's model based values (excluding
        %NaNs)
        diff_evidence_streams = mean(MB_same_nback','omitnan')' - mean(MB_diff_nback','omitnan')';
        
        
        %insert descriptive variables to enable model fitting
        trials = [1:n_trials]';
        this_sub = repmat(subs, 1, n_trials)';
        rewards_s1 = rewards(:, 1); %1/-1
        rewards_s2 = rewards(:, 2);% 0
        
        %save simulated behaviour in a table; to be used for model fitting
        dat_tab = table(this_sub, trials, correlated_block, gem_presented, market_presented, probs, rewards_s1, rewards_s2,...
            prob_choose_left,choose_left, state, outcome, ...
            MB_same_1_back, MB_same_2_back, MB_same_3_back, MB_same_4_back, MB_same_5_back,...
            MB_diff_1_back, MB_diff_2_back, MB_diff_3_back, MB_diff_4_back, MB_diff_5_back, diff_evidence_streams, forced_trial);
        
        if subs==1
            dat_tab_all = dat_tab;
        else
            dat_tab_all = [dat_tab_all; dat_tab]; %concatenate subject files 
        end
        
        
    end
    
    %% model fitting to the analysed simulated data
    dat_tab_all.correlated_block(dat_tab_all.correlated_block == 2) = -1; 
    
    %model formula for lagged regression
    m3 = fitglme(dat_tab_all((dat_tab_all.forced_trial==0), :), 'choose_left ~ diff_evidence_streams*correlated_block + (1 + diff_evidence_streams*correlated_block |this_sub)', 'Distribution', 'binomial'); 

    stats_difference.stat_estimate(sims, :) = m3.Coefficients.Estimate(2:end)';
    stats_difference.stat_comp_tstat(sims, :) = m3.Coefficients.tStat(2:end)';
    stats_difference.stat_comp_pval(sims, :) = m3.Coefficients.pValue(2:end)';
    stats_difference.stat_SE(sims, :) = m3.Coefficients.SE(2:end);
    
    stats_difference.estimate_interaction(sims, 1) = m3.Coefficients.Estimate(4);
    stats_difference.tstat_interaction(sims, 1) = m3.Coefficients.tStat(4);

end

hist(stats_difference.tstat_interaction);
hist(stats_difference.estimate_interaction);

quantile(stats_difference.tstat_interaction, [0.025, 0.0975])
quantile(stats_difference.estimate_interaction, [0.025, 0.0975])

csvwrite('tstat_permutation.csv', stats_difference.tstat_interaction)
csvwrite('estimate_permutation.csv', stats_difference.estimate_interaction)


