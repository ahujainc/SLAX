%-------------------------------------------------------------------------
% Function Name: fun_grid_gen
% Date last modified: December 12, 2018
% Author: Vishal Ahuja

% PURPOSE: This function takes seat of grid points (states1) and simplices
% (Del_T1) and returns the refined set of states (Del_T2) and simplices
% (Del_T2) - essentially increasing the grid resoution (Z) by 1. The
% refinement is done such that the center points of each of the existing
% simplices are added to the grid state space.
%-------------------------------------------------------------------------
function [states2 Del_T2] = fun_grid_gen(states1,Del_T1)
global d
states2=states1;
Del_T2=Del_T1;
ix=size(states1,1);
for i=1:length(Del_T1)
    simplex_states=[];
    for j=1:d+1
        simplex_states=[simplex_states; states1(Del_T1(i,j),:)];
    end
    temp_state = mean(simplex_states);
    states2=[states2;temp_state];
    ix=ix+1;
    curr_simplex = Del_T1(i,:);
    add_simplex = [];
    for k=1:d+1
        temp_simplex = curr_simplex;
        ix_val = curr_simplex(k);
        temp_simplex(temp_simplex==ix_val) = ix;
        add_simplex = [add_simplex; temp_simplex];
    end
    Del_T2(1,:) = [];
    Del_T2 = [Del_T2; add_simplex];
end