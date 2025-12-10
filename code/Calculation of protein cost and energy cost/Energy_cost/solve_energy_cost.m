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
ub_of_compositions = [0,1000,100,150,1000,1000,1000,...
    1000,200,1000,130,1000,0,130,...
    100,400,400,200,1000,100,1000,...
    90,1000,170,1000,1000,300,170,...
    30,100,170,1000,1000]'; 

uptake_flux_id_table = struct2table(rxns_culture_media);
uptake_flux_ids = table2array(uptake_flux_id_table);
bux(uptake_flux_ids) = ub_of_compositions;

%% Maximize the biomass production for 'reference state' 
c_biomass = ecModel1.c(1:n_true_rxns);
% blc = zeros(n_true_mets,1);
% buc = zeros(n_true_mets,1);
% 
% % [18354-ala, 18361-arg, 18235-asn, 18236-asp, c-cys, 18358-gln, 8850-glu, 18351-gly, 
% % 18273-pro, 18275-ser, 18360-tyr]
% bux_ban_aa_of_interest = bux;
% bux_ban_aa_of_interest(18360) = 0; % (*******)
% 
% [r, res] = mosekopt('maximize echo(2)', struct('c', c_biomass, 'a', ...
%     S_rm_prot_pool, 'blc', blc, 'buc', buc, 'blx', blx, 'bux', bux_ban_aa_of_interest));
% biomass_optval = res.sol.bas.pobjval;
% disp(biomass_optval);

%% Maximize the biomass production for 'case state'
% S_rm_prot_pool_ds = S_rm_prot_pool;
% ds = -0.1;
% % [552-ala, 658-arg, 553-asn, 149-asp, 555-cys, 559-gln,123-glu, 561-gly, 
% % 2493-pro, 563-ser, 2091-tyr]
% S_rm_prot_pool_ds(2493, biomass_id) = S_rm_prot_pool_ds(2493, biomass_id) + ds; % (*******)
% 
% [r1, res1] = mosekopt('maximize echo(2)', struct('c', c_biomass, 'a', ...
%     S_rm_prot_pool_ds, 'blc', blc, 'buc', buc, 'blx', blx, 'bux', bux_ban_aa_of_interest));
% biomass_optval1 = res1.sol.bas.pobjval;

%% Minimize the energy cost with fixed biomass flux
ATP_indices_in_metNames = find(strcmp(ecModel1.metNames, 'Adenosine Triphosphate')); 
ATP_merged_row_of_S = sum(S_rm_prot_pool(ATP_indices_in_metNames,:), 1);
c = ATP_merged_row_of_S';
c(c > 0) = 0;
c = c * -1;

blc = zeros(n_true_mets,1);
buc = zeros(n_true_mets,1);

blx_add_biomass_lb = blx;
biomass_id = find(c_biomass>0);
blx_add_biomass_lb(biomass_id) = 0.7 * 20.1431;

% [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% 18273-pro, 18275-ser, 18360-tyr]
bux_ban_aa_of_interest = bux;
bux_ban_aa_of_interest(18360) = 0; % (*******)

[r_ref, res_ref] = mosekopt('minimize echo(2)', struct('c', c, 'a', S_rm_prot_pool,...
    'blc', blc, 'buc', buc, 'blx', blx_add_biomass_lb, 'bux', bux_ban_aa_of_interest));

result_ref = res_ref.sol.bas.pobjval;
v_ref = res_ref.sol.bas.xx;
disp(result_ref);

%% Repeat the steps above when a slight ds is added in the stoichiometric coefficient of
% the AA of interest in the biomass composition 
S_rm_prot_pool_ds = S_rm_prot_pool;
ds = -0.01;
% [552-ala, 658-arg, 553-asn, 149-asp, 555-cys, 559-gln,123-glu, 561-gly, 
% 2493-pro, 563-ser, 2091-tyr]
S_rm_prot_pool_ds(2091, biomass_id) = S_rm_prot_pool_ds(2091, biomass_id) + ds; % (*******)

[r_case, res_case] = mosekopt('minimize echo(2)', struct('c', c, 'a', S_rm_prot_pool_ds, ...
    'blc', blc, 'buc', buc, 'blx', blx_add_biomass_lb, 'bux', bux_ban_aa_of_interest));

result_case = res_case.sol.bas.pobjval;
v_case = res_case.sol.bas.xx;
disp(result_case);

cost_approx = (result_case - result_ref)/abs(ds);
disp(cost_approx);

%% find sources (case - reference)
indices_targeted_metabolite = find(strcmp(ecModel1.metNames, 'L-Tyrosine')); % (*******)
row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
v_difference = v_case - v_ref;
source_metabolite = row_targeted_metabolite' .* v_difference;
