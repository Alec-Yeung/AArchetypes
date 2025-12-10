%% Human
w = 5.77;                 % hsa: 5.77、sce: 0.26、eco: 0.54
d_vec = linspace(0, 6, 10000);  
nd = numel(d_vec);
c_ext = 1;

a_opt = zeros(nd,1);
A_opt = zeros(nd,1);
obj_opt = zeros(nd,1);
ext_cost = zeros(nd,1);    

lb = [0; 0];             
ub = [1; Inf];           
opts = optimoptions('fmincon','Display','off','Algorithm','sqp-legacy');

for k = 1:nd
    d = d_vec(k);
    fun = @(x) benefit_Human(x, d);           
    nonlcon = @(x) budget_constraint_Human(x, w);


    x0 = [0.5; 0.5];                    
    [xbest, fval] = fmincon(fun, x0, [],[],[],[], lb, ub, nonlcon, opts);

    a_opt(k) = xbest(1);
    A_opt(k) = xbest(2);
    obj_opt(k) = -fval;                 

    ext_cost(k) = external_cost_Human(a_opt(k), c_ext);
end

% 图1：a ~ d
figure;
plot(d_vec, a_opt, 'LineWidth', 1);
xlabel('d');
ylabel('optimal a');
title('Human');
grid off;
ylim([0 1]);   

% 图2：external cost ~ d
figure;
plot(d_vec, ext_cost, 'LineWidth', 1);
xlabel('d');
ylabel('external cost');
title('Human');
grid off;
ylim([0 1]);   


%% 分别求两个w（hsa和sce）下的优化问题，把两条ec~d曲线画在一张图上
w_list = [5.77, 0.26];              
colors = lines(numel(w_list));     
d_vec = linspace(0, 6, 10000);
nd = numel(d_vec);
c_ext = 1;

ext_cost_all = zeros(nd, numel(w_list));

for i = 1:numel(w_list)
    w = w_list(i);

    a_opt = zeros(nd,1);
    A_opt = zeros(nd,1);
    obj_opt = zeros(nd,1);
    ext_cost = zeros(nd,1);

    lb = [0; 0];
    ub = [1; Inf];
    opts = optimoptions('fmincon','Display','off','Algorithm','sqp-legacy');

    for k = 1:nd
        d = d_vec(k);
        fun = @(x) benefit_Human(x, d);
        nonlcon = @(x) budget_constraint_Human(x, w);
        x0 = [0.5; 0.5];
        [xbest, fval] = fmincon(fun, x0, [],[],[],[], lb, ub, nonlcon, opts);

        a_opt(k) = xbest(1);
        A_opt(k) = xbest(2);
        obj_opt(k) = -fval;
        ext_cost(k) = external_cost_Human(a_opt(k), c_ext);
    end

    ext_cost_all(:, i) = ext_cost;
end

figure;
hold on;
for i = 1:numel(w_list)
    plot(d_vec, ext_cost_all(:, i), 'LineWidth', 1.5, 'Color', colors(i,:));
end
xlabel('d');
ylabel('external cost');
title('Human');
legend(arrayfun(@(w) sprintf('w = %.2f', w), w_list, 'UniformOutput', false), ...
       'Location', 'best');
grid off;
ylim([0 1]);
box on
hold off;

%% 分别求两个w（hsa和sce）下的优化问题，把两条a~d曲线画在一张图上
w_list = [5.77, 0.26];              
colors = lines(numel(w_list));     
d_vec = linspace(0, 6, 10000);
nd = numel(d_vec);
c_ext = 1;

a_all = zeros(nd, numel(w_list));

for i = 1:numel(w_list)
    w = w_list(i);

    a_opt = zeros(nd,1);
    A_opt = zeros(nd,1);
    obj_opt = zeros(nd,1);
    ext_cost = zeros(nd,1);

    lb = [0; 0];
    ub = [1; Inf];
    opts = optimoptions('fmincon','Display','off','Algorithm','sqp-legacy');

    for k = 1:nd
        d = d_vec(k);
        fun = @(x) benefit_Human(x, d);
        nonlcon = @(x) budget_constraint_Human(x, w);
        x0 = [0.5; 0.5];
        [xbest, fval] = fmincon(fun, x0, [],[],[],[], lb, ub, nonlcon, opts);

        a_opt(k) = xbest(1);
        A_opt(k) = xbest(2);
        obj_opt(k) = -fval;
        ext_cost(k) = external_cost_Human(a_opt(k), c_ext);
    end

    a_all(:, i) = a_opt;
end

% 绘图
figure;
hold on;
for i = 1:numel(w_list)
    plot(d_vec, a_all(:, i), 'LineWidth', 1.5, 'Color', colors(i,:));
end
xlabel('d');
ylabel('a');
title('Human');
legend(arrayfun(@(w) sprintf('w = %.2f', w), w_list, 'UniformOutput', false), ...
       'Location', 'best');
grid off;
ylim([0 1]);
box on
hold off;


%%
function p = benefit_Human(x,d)
a = x(1);
A = x(2);
pA = LLPS_prob_Human(a);

benefit = A*(1+(d)*pA);
p = -benefit;
end

function p = benefit_Yeast(x,d)
a = x(1);
A = x(2);
pA = LLPS_prob_Yeast(a);

benefit = A*(1+(d)*pA);
p = -benefit;
end

function p = LLPS_prob_Human(x)
p = 1 / (1 + exp(6.3296-10.5557*x));
end

function p = LLPS_prob_Yeast(x)
p = 1 / (1 + exp(7.0094-11.5778*x));
end

function [c,ceq] = budget_constraint_Human(x,w)
a = x(1);
A = x(2);

c = (w*a+(1-a))*A - 1;
ceq = [];
end

function [c,ceq] = budget_constraint_Yeast(x,w)
a = x(1);
A = x(2);

c = (w*a+(1-a))*A - 1;
ceq = [];
end

function cost = external_cost_Human(a, c)
pA = LLPS_prob_Human(a);
cost = c * pA;
end

function cost = external_cost_Yeast(a, c)
pA = LLPS_prob_Yeast(a);
cost = c * pA;
end




