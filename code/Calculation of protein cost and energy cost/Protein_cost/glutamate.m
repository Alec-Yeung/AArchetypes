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
ub_of_compositions = [102,56,14,14,1000,1000,1000,...
    1000,15,22,22,1000,0,127,...
    22,52,88,63,1000,20,1000,...
    44,1000,100,1000,1000,75,51,...
    10,29,75,1000,1000]'; 

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

%% Minimize the protein cost with fixed biomass flux
w_modify = w;
massive_w_ids = find(w > 1e+4);
w_modify(massive_w_ids) = 1e+4;
c_protein_cost = w_modify;

blc = zeros(n_true_mets,1);
buc = zeros(n_true_mets,1);

blx_add_biomass_lb = blx;
biomass_id = find(c_biomass>0);
blx_add_biomass_lb(biomass_id) = 5;

% [18354-ala, 18361-arg, 18235-asn, 18236-asp, 18353-cys, 18358-gln, 8850-glu, 18351-gly, 
% 18273-pro, 18275-ser, 18360-tyr]
bux_ban_aa_of_interest = bux;
% bux_ban_aa_of_interest(8850) = 0; % (*******)

[r_ref, res_ref] = mosekopt('minimize echo(2)', struct('c', c_protein_cost, 'a', ...
    S_rm_prot_pool, 'blc', blc, 'buc', buc, 'blx', blx_add_biomass_lb, 'bux', ...
    bux_ban_aa_of_interest));
result_ref = res_ref.sol.bas.pobjval;
v_ref = res_ref.sol.bas.xx;
disp(result_ref);

%% Repeat the steps above when a slight ds is added in the stoichiometric coefficient of 
% the AA of interest in the biomass composition 
S_rm_prot_pool_ds = S_rm_prot_pool;
ds = -0.1;
% [552-ala, 658-arg, 553-asn, 149-asp, 555-cys, 559-gln,123-glu, 561-gly, 
% 2493-pro, 563-ser, 2091-tyr]
S_rm_prot_pool_ds(123, biomass_id) = S_rm_prot_pool_ds(123, biomass_id) + ds; % (*******)

indices_targeted_metabolite = find(strcmp(ecModel1.metNames, 'L-Glutamate')); 
row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
rxns_glu_engaged = find(row_targeted_metabolite ~= 0)';
targeted_rxn = (rxns_glu_engaged == 14448); 
rxns_glu_engaged_that_remain_unchanged_flux = rxns_glu_engaged(~targeted_rxn);
blx_case = blx_add_biomass_lb;
blx_case(rxns_glu_engaged_that_remain_unchanged_flux) = v_ref(rxns_glu_engaged_that_remain_unchanged_flux); 
bux_case = bux_ban_aa_of_interest;
bux_case(rxns_glu_engaged_that_remain_unchanged_flux) = v_ref(rxns_glu_engaged_that_remain_unchanged_flux); 
% case状态下，glu参与的反应的flux (除14448外) 与ref状态保持一致

[r_case, res_case] = mosekopt('minimize echo(2)', struct('c', c_protein_cost, 'a', ...
    S_rm_prot_pool_ds, 'blc', blc, 'buc', buc, 'blx', blx_case, 'bux', bux_case));

result_case = res_case.sol.bas.pobjval;
v_case = res_case.sol.bas.xx;
disp(result_case);

cost_approx = (result_case - result_ref)/abs(ds);
disp(cost_approx);

%% find sources (case - reference)
indices_targeted_metabolite = find(strcmp(ecModel1.metNames, 'L-Glutamate')); % (*******)
row_targeted_metabolite = sum(S_rm_prot_pool(indices_targeted_metabolite,:), 1);
v_difference = v_case - v_ref;
source_metabolite = row_targeted_metabolite' .* v_difference;
distribution_of_engergy_cost = c_protein_cost .* v_difference;

%% 找出w_modify和v_difference均不为0的反应

% 找出 w 和 v 同时不为 0 的索引
rxns_contribute_to_cost = find(abs(w_modify) > 1e-4 & abs(v_difference) > 1e-4);

% 获取对应的 w 和 v 元素，并转换为全矩阵
w_nonzero = full(w_modify(rxns_contribute_to_cost));
v_nonzero = full(v_difference(rxns_contribute_to_cost));

% 使用 arrayfun 创建字符串乘积
char_products = arrayfun(@(x, y) sprintf('%g*%g', x, y), w_nonzero, v_nonzero, 'UniformOutput', false);

% 创建 Cell Array
result = [num2cell(rxns_contribute_to_cost), char_products];
