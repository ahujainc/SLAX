# Matlab Code for implementing SLAX
# An Approximation Approach for Response Adaptive Clinical Trial Design (https://ssrn.com/abstract=3212148)
by Vishal Ahuja and John R. Birge

Main function file: "mainfun_SLAX" - this function is the main function and calls the following functions:
1) fun_ap_sim_B
2) fun_ap_imp_sim_B

Notes:
- Both the functions "fn_ap_sim_B" and "fun_ap_imp_sim_B" call the function "fun_VFap_sim"
- If needed, the function "fun_grid_gen" could be called from the "mainfun_SLAX" function (for further grid refinement). Additional, details are provided in each of the functions.
