%-------------------------------------------------------------------------
% Function Name: fun_ap_sim_B
% Date last modified: December 12, 2018
% Author: Vishal Ahuja
% PURPOSE: Compute and store approximate values to go for each grid point
% This function calls the following functions: fun_VFap_sim 
%-------------------------------------------------------------------------

function VF=fun_ap_sim_B
global VF_sim states num_states ab pat T I actions Del_T alpha srun2 B numgen H

% ix = T+1; % set the time index; initialze at 1 (we dont need to evaluate
% for all time periods)

ix=H+2;  % set the time index; initialze at 1
VF_sim = zeros(num_states,ix); % Initialize the value function
state_init = zeros(1,size(states,2)); % A single initial state: [0 0 0]

func_txt = str2func('fun_VFap_sim'); % string to function converter

% Loop through all time periods (ignore t=H since no decision needs to be made)
for t=ix-1:-1:1 
    % set the states
    if t>1
        n_Total = num_states;
        curr_states = states;
    else
        n_Total = 1; %if t=0, then there is only a single state: [0 0 0]
        curr_states = state_init;
    end
    % loop through all the states
    for n=1:n_Total        
        
        curr_state = curr_states(n,:);  % current state being evaluated
        nt=(t-1)*pat;                   % total number of patients
        
        n1=curr_state(1)*nt; % number of patients allocated to A so far
        n2=nt-n1;            % number of patients allocated to B so far
        s1=n1*curr_state(2); % number of successes on A so far
        s2=n2*curr_state(3); % number of successes on B so far
        f1=n1-s1;            % number of failures on A so far
        f2=n2-s2;            % number of failures on B so far
        p1 = (ab(1,1)+s1)/(sum(ab(1,:))+n1); % Expected Prob of success for A
        p2 = (ab(2,1)+s2)/(sum(ab(2,:))+n2); % Expected Prob of success for B
       
        
        if t==ix-1
            % if last decision period, pick the Tx with the highest success
            % probability (greedy) and evaluate VF accordingly
            VF_sim(n,t)=pat*max(p1,p2);
        else
            % simulate to get value function
            Esuc_act = zeros(actions,1);  % initialize mean number of successes
            Stdev_act = zeros(actions,1); % initialize stdev of number of successes
            UCB = zeros(actions,1); % weighted objective function
            
            
            for a=1:B:actions  % loop through actions in blocks
                Esuc=NaN(srun2,1); % initialize mean number of successes (inner loop)
                n1_act = a-1;
                n2_act = pat-n1_act;
                % Use simulation to generate optimal action
                for r2=1:srun2
                    s1_act = binornd(n1_act,p1); % successes on Tx A
                    s2_act = binornd(n2_act,p2); % successes on Tx B
                    f1_act = n1_act-s1_act;      % failures on Tx A
                    f2_act = n2_act-s2_act;      % failures on Tx B
                    
                    nextst = [s1 f1 s2 f2] + [s1_act f1_act s2_act f2_act]; 
                    
                    % V(s,a)=sum of successes in this period + value-to-go for next state
                    % which is interpolated (if a non-grid point) using the function "fun_VFap_sim"
                    Esuc(r2) = s1_act + s2_act + func_txt(nextst,t);
                end
                Esuc_act(a)= mean(Esuc);
                Stdev_act(a)= std(Esuc);
                
                % Objective function is the weighted sum of mean and stddev
                UCB(a)=(1-alpha)*Esuc_act(a)+alpha*Stdev_act(a);
            end
            
            opt_act=find(UCB==max(UCB),1,'first'); % optimal action is the one with highest obj. 
            VF_sim(n,t) = Esuc_act(opt_act); % VF of grid state n at time t
        end
    end
end
VF=VF_sim;
end