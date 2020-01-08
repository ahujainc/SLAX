%-------------------------------------------------------------------------
% Function Name: fun_ap_imp_sim_B
% Date last modified: December 12, 2018
% Author: Vishal Ahuja
% PURPOSE: Implements SLAX on all the possible states
% This function calls the following functions: fun_VFap_sim 
%-------------------------------------------------------------------------
function mixed_opt=fun_ap_imp_sim_B
global VF_sim ab pat T I actions Del_T alpha srun1 srun2 numgen B H

ix = T+1;  % creates the index

E_num_suc=zeros(srun1,1);% initialize the number of successes (outer loop)
func_txt = str2func('fun_VFap_sim'); % calls the function for bounds

for r1=1:srun1 % loop thr
    for t=1:ix
        
        % initialize the current state
        if t==1
            curr_state = [0 0 0 0];
            val2go = 0;
        else
            curr_state = curr_state + [s1_opt f1_opt s2_opt f2_opt];
        end
        s1 = curr_state(1); f1 = curr_state(2); n1 = s1+f1;
        s2 = curr_state(3); f2 = curr_state(4); n2 = s2+f2;
        p1 = (ab(1,1)+s1)/(sum(ab(1,:))+n1);
        p2 = (ab(2,1)+s2)/(sum(ab(2,:))+n2);
        
        if t<ix-1 % CHANGE THIS TO ADJUST FOR H
            % see file "fun_ap_sim_B" for details on the variables below
            if t<H+1
                Esuc_act = zeros(actions,1);
                Stdev_act = Esuc_act;
                UCB = Esuc_act;
                for a=1:B:actions % loop through actions in blocks
                    Esuc=NaN(srun2,1); % this is the outer loop
                    n1_act = a-1;
                    n2_act = pat-n1_act;
                    % Use simulation to generate optimal action
                    for r2=1:srun2
                        s1_act = binornd(n1_act,p1);
                        s2_act = binornd(n2_act,p2);
                        f1_act = n1_act-s1_act;
                        f2_act = n2_act-s2_act;
                        nextst = [s1 f1 s2 f2] + [s1_act f1_act s2_act f2_act];
                        % Calls the function "fun_VFap_sim"
                        Esuc(r2) = s1_act + s2_act + func_txt(nextst,t);
                    end
                    Esuc_act(a)= mean(Esuc);
                    Stdev_act(a)= std(Esuc);
                    UCB(a)=(1-alpha)*Esuc_act(a)+alpha*Stdev_act(a);
                end
                opt_act=find(UCB==max(UCB),1,'first');
                n1_opt = opt_act-1;
                % for the remainder (T-H) periods use Greedy heuristic
            else
                if p1>=p2
                    n1_opt=pat;
                else
                    n1_opt=0;
                end
            end
            n2_opt = pat-n1_opt;
      
            % simulate number of successes on each arm
            s1_opt = binornd(n1_opt,p1);
            s2_opt = binornd(n2_opt,p2);            
            f1_opt = n1_opt-s1_opt;
            f2_opt = n2_opt-s2_opt;
            val2go = val2go+s1_opt+s2_opt;
            
            % use Greedy Heuristics for the terminal decision stage
        elseif t==ix-1
            val2go = val2go+pat*max(p1,p2);
        elseif t==ix
            % This is the value function for a given (outer loop) run
            E_num_suc(r1) = val2go;
        end
    end
end
mixed_opt=[mean(E_num_suc), std(E_num_suc)/sqrt(srun1)];
end
