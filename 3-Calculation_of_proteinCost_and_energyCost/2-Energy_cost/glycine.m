%% Find reactions that produces metabolites from nothing and block fluxes through them
load('w.mat');
load('rxns_culture_media.mat');
load('ecModel.mat');

n_true_mets = 8399;
n_true_rxns = 20222;
S_rm_prot_pool = ecModel1.S(1:n_true_mets, 1:n_true_rxns);
creater_rxn_ids = find(min(S_rm_prot_pool)==0);
blx = zeros(n_true_rxns,1);
bux = 1000*ones(n_true_rxns,1);
bux(creater_rxn_ids) = 0;

%% Allow fluxes through uptake reactions
ub_of_compositions = [97,53,13,13,1000,1000,1000,...
    1000,13,22,22,1000,0,129,...
    25,51,91,72,1000,23,1000,...
    42,1000,105,1000,1000,70,50,...
    10,31,69,1000,1000]'; 

uptake_flux_id_table = struct2table(rxns_culture_media);
uptake_flux_ids = table2array(uptake_flux_id_table);
bux(uptake_flux_ids) = ub_of_compositions;

%% Maximize the biomass production for 'reference state' 
% c_biomass = ecModel1.c(1:n_true_rxns);
% 
% blc = zeros(n_true_mets,1);
% buc = zeros(n_true_mets,1);
% 
% blx0 = blx;
% bux0 = bux;
% % [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% % 18273-pro, 18275-ser, 18360-tyr]
% bux0(18273) = 0; % (*******)
% 
% [r0, res0] = mosekopt('maximize echo(2)', struct('c', c_biomass, 'a', S_rm_prot_pool, ...
%     'blc', blc, 'buc', buc, 'blx', blx0, 'bux', bux0));
% biomass_optval_ref = res0.sol.bas.pobjval;
% v0 = res0.sol.bas.xx;
% disp(biomass_optval_ref);
% 
% % biomass_ref_Pro = 20.4427

%% Minimize the energy cost with fixed biomass flux-reference state 
ATP_indices_in_metNames = find(strcmp(ecModel1.metNames, 'Adenosine Triphosphate')); 
ATP_merged_row_of_S = sum(S_rm_prot_pool(ATP_indices_in_metNames,:), 1);
c = ATP_merged_row_of_S';
c(c > 0) = 0;
c = c * -1;
% w_modify = w;
% massive_w_ids = find(w > 1e+4);
% w_modify(massive_w_ids) = 1e+4;
% c_protein_cost = w_modify;

blc = zeros(n_true_mets,1);
buc = zeros(n_true_mets,1);

blx_ref = blx;
bux_ref = bux;
c_biomass = ecModel1.c(1:n_true_rxns);
biomass_id = find(c_biomass>0);
blx_ref(biomass_id) = 5; % (*******)
% [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% 18273-pro, 18275-ser, 18360-tyr]
bux_ref(18351) = 0; % (*******)

bux_ref(1810) = 500;

[r_ref, res_ref] = mosekopt('minimize echo(2)', struct('c', c, 'a', S_rm_prot_pool,...
    'blc', blc, 'buc', buc, 'blx', blx_ref, 'bux', bux_ref));

result_ref = res_ref.sol.bas.pobjval;
v_ref = res_ref.sol.bas.xx;
disp(result_ref);

%% Maximize the biomass production for 'case state'
% c_biomass = ecModel1.c(1:n_true_rxns);
% 
% S_rm_prot_pool_ds = S_rm_prot_pool;
% ds = -0.1;
% biomass_id = find(c_biomass>0);
% % [552-ala, 658-arg, 553-asn, 149-asp, 555-cys, 559-gln,123-glu, 561-gly, 
% % 2493-pro, 563-ser, 2091-tyr]
% S_rm_prot_pool_ds(561, biomass_id) = S_rm_prot_pool_ds(561, biomass_id) + ds; % (*******)
% 
% blc = zeros(n_true_mets,1);
% buc = zeros(n_true_mets,1);
% 
% blx1 = blx;
% bux1 = bux;
% % [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% % 18273-pro, 18275-ser, 18360-tyr]
% bux1(18351) = 0; % (*******)
% 
% indices_targeted_metabolite = strcmp(ecModel1.metNames, 'Glycine'); % (*******)
% row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
% rxns_AA_engaged = find(row_targeted_metabolite ~= 0)';
% targeted_rxn = (rxns_AA_engaged == 1810); % (*******)
% rxns_AA_engaged_that_remain_unchanged_flux = rxns_AA_engaged(~targeted_rxn);
% blx1(rxns_AA_engaged_that_remain_unchanged_flux) = v_ref(rxns_AA_engaged_that_remain_unchanged_flux); 
% bux1(rxns_AA_engaged_that_remain_unchanged_flux) = v_ref(rxns_AA_engaged_that_remain_unchanged_flux); 
% 
% [r1, res1] = mosekopt('maximize echo(2)', struct('c', c_biomass, 'a', S_rm_prot_pool_ds, ...
%     'blc', blc, 'buc', buc, 'blx', blx1, 'bux', bux1));
% 
% biomass_optval_case = res1.sol.bas.pobjval;
% v1 = res1.sol.bas.xx;
% disp(biomass_optval_case)

%% Minimize the energy cost with fixed biomass flux-case state 
S_rm_prot_pool_ds = S_rm_prot_pool;
ds = -0.1;
% [552-ala, 658-arg, 553-asn, 149-asp, 555-cys, 559-gln,123-glu, 561-gly, 
% 2493-pro, 563-ser, 2091-tyr]
S_rm_prot_pool_ds(561, biomass_id) = S_rm_prot_pool_ds(561, biomass_id) + ds; % (*******)

blx_case = blx;
bux_case = bux;
% [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% 18273-pro, 18275-ser, 18360-tyr]
bux_case(18351) = 0; % (*******)

indices_targeted_metabolite = strcmp(ecModel1.metNames, 'Glycine'); % (*******)
row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
rxns_AA_engaged = find(row_targeted_metabolite ~= 0)';
targeted_rxn = (rxns_AA_engaged == 1810); % (*******)
rxns_AA_engaged_that_remain_unchanged_flux = rxns_AA_engaged(~targeted_rxn);

blx_case(rxns_AA_engaged_that_remain_unchanged_flux) = v_ref(rxns_AA_engaged_that_remain_unchanged_flux); 
bux_case(rxns_AA_engaged_that_remain_unchanged_flux) = v_ref(rxns_AA_engaged_that_remain_unchanged_flux); 
% case状态下，gly参与的反应的flux (除1810外) 与ref状态保持一致 % (*******)

[r_case, res_case] = mosekopt('minimize echo(2)', struct('c', c, 'a', S_rm_prot_pool_ds,...
    'blc', blc, 'buc', buc, 'blx', blx_case, 'bux', bux_case));

result_case = res_case.sol.bas.pobjval;
v_case = res_case.sol.bas.xx;
disp(result_case);

cost_approx = (result_case - result_ref)/abs(ds);
disp(cost_approx);

%% find sources (case - reference)
indices_targeted_metabolite = find(strcmp(ecModel1.metNames, 'Glycine')); % (*******)
row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
v_difference = v_case - v_ref;
source_metabolite = row_targeted_metabolite' .* v_difference;
distribution_of_engergy_cost = c .* v_difference;

%% 
flux_difference_of_rxns_AA_engaged = zeros(length(rxns_AA_engaged), 2);
flux_difference_of_rxns_AA_engaged(:, 1) = rxns_AA_engaged;
flux_difference_of_rxns_AA_engaged(:, 2) = v_difference(rxns_AA_engaged);
filter = flux_difference_of_rxns_AA_engaged(:, 2) > 0.00001 | ...
    flux_difference_of_rxns_AA_engaged(:, 2) < -0.00001;
flux_difference_of_rxns_AA_engaged = flux_difference_of_rxns_AA_engaged(filter, :);
