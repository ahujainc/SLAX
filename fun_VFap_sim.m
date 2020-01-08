%-------------------------------------------------------------------------
% Function Name: fun_VFap_sim
% Date last modified: November 12, 2018
% Author: Vishal Ahuja
% This function interpolates the value function of a non-grid state using
% Kuhn's triangulation as laid out in Munos and Moore (1999). Further, the
% val function = min(max(interpolated value, lower bound), upper % bound))
%-------------------------------------------------------------------------
function val2go=fun_VFap_sim(svals,t)
global VF_sim states num_states ab pat T I Del_T ub_list

% we first establish the state whose value is to be interpolated in terms
s1=svals(1); % successes on Tx A
f1=svals(2); % successes on Tx B
s2=svals(3); % failures on Tx A
f2=svals(4); % failures on Tx B
n1=s1+f1;    % Total patients on Tx A
n2=s2+f2;    % Total patients on Tx B
ab_post = ab + [s1 f1; s2 f2]; % (alpa,beta) parameters 

% ----------STEP 1: Interpolation-------
% Establish the state based on success/failures in each Tx

fr_1=n1/(n1+n2); fr_1s=s1/n1; fr_2s=s2/n2;
% if the denominator is zero, then the fraction is zero
if(isnan(fr_1s) && (s1==0) && (n1==0))
    fr_1s=0;
end
if(isnan(fr_2s) && (s2==0) && (n2==0)) 
    fr_2s=0; 
end

pquery = [fr_1 fr_1s fr_2s]; % state that needs to be interpolated
points = states; % all the available grid states

% The following set of code finds the simplex that contains the nongrid
% point and solves the linear system of eq-ns to find barycentric weights
for i=1:size(Del_T,1)
    psimplex=points(Del_T(i,:),:);
    A = [psimplex'; 1 1 1 1];
    B = [pquery';1];
    barcoords = A\B; % solve a system of linear equations
    if all(barcoords>=-10^-6)
        interp_simplex =i;
        interp_barcoords=barcoords';
        break;
    end
end
f_val_ip=interp_barcoords*VF_sim(Del_T(interp_simplex,:)',t+1);

% ----------STEP 2: Lower bounds = max (Ep1,Ep2) (greedy algorithm)-----
Ep1=ab_post(1,1)/sum(ab_post(1,:));
Ep2=ab_post(2,1)/sum(ab_post(2,:));
f_val_lb=max(Ep1,Ep2)*(T-t)*pat;

% ----------STEP 3: Upper bounds:max(p1,p2)=Emax{pr(p1>p2),Pr(p2>p1)}------ 
ab2=[ab_post(1,:) ab_post(2,:)];
% first check if the file ub_list is uploaded
if isempty(ub_list)==1
    indx=[];
else
    indx=find(ub_list(:,1)==ab2(1) & ub_list(:,2)==ab2(2) & ub_list(:,3)==ab2(3) & ub_list(:,4)==ab2(4),1,'last');
end
% [~,indx]=ismember(ab2,ub_list(:,1:4),'rows');
if isempty(indx)==1
    prob_ub1=integral(@(x)x.*betapdf(x,ab_post(2,1),ab_post(2,2)).*betacdf(x,ab_post(1,1),ab_post(1,2)),0,1); %E[p1|p1>p2]
    prob_ub2=integral(@(y)y.*betapdf(y,ab_post(1,1),ab_post(1,2)).*betacdf(y,ab_post(2,1),ab_post(2,2)),0,1); %E[p2|p2>p1]
    f_val_ub=(prob_ub1+prob_ub2)*(T-t)*pat;
else
    prob_ub=ub_list(indx,5);
    f_val_ub=prob_ub*(T-t)*pat;
end

% prob_ub=ub_list(indx,5);
% f_val_ub=prob_ub*(T-t)*pat;

% VF = min(max(interpolated value, lower bound), upper % bound))
% ----------------------------------------------------------------
f_val = min(max(f_val_ip,f_val_lb),f_val_ub);
val2go=f_val;
