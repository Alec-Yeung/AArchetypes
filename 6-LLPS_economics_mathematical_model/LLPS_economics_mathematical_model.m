% Mathematical model for amino acid composition optimization
% ------------ Basic setting of the model --------------- %
% Two amino acids, x and y, are combined into a protein P with abundance A
% ax + (1-a)y -> P
% The amino acid x increases the chance of protein LLPS, therefore
% increasing the functional payoff by dp(a) for one unit of protein
% p(a) was estimated based on protein LLPS database: p(a) =
% 1/(1+exp(6.3-10.6a))
% Cost of x is w, cost of y is 1
% w is a species-specific LLPS cost ratio
% Total budget (w*a+(1-a)*A cannot exceed 1
% The functional payoff function P_func(a,A)=(1+dp(a))A
% Information capacity term P_info(a)=-alog2(a/7)-(1-a)log2((1-a)/13)
% Overall payoff to be optimized: P(a,A) = P_func(a,A)+beta*P_info(a)

%% Plot w for H.sapiens (human) and S.cerevisiae (yeast)
w_human = 5.77;
w_yeast = 0.26;

figure;
bar([1 2],[w_human w_yeast],'FaceColor',[233 113 50]/256,'EdgeColor',[1 1 1]);
ylabel("w");
xticklabels(["H. sapiens" "S. cerevisiae"]);xtickangle(45);
title("LLPS cost ratio (w)");
hold on;
plot([0 3],[1 1],':','Color',[0.5 0.5 0.5]);
xlim([0.2 2.8]);

%% Read and visualize real proteomics data
a_data_table = readtable("Ratio of IDR-associated AAs.xlsx");
a_data_human = a_data_table{ismember(a_data_table.group,'Homo sapiens'),"value"};
a_data_yeast = a_data_table{ismember(a_data_table.group,'Saccharomyces cerevisiae'),"value"};
figure;
subplot(2,1,1);
histogram(a_data_human,0.2:0.01:0.7, ...
    'FaceColor',[233 113 50]/256,'EdgeColor',[1 1 1]);
xlim([0.2 0.7]);title("Distribution of a, H. sapiens",'FontWeight','Normal');
xlabel("a");ylabel("#Proteins");
ax = gca;
ax.LineWidth = 1;
subplot(2,1,2);
histogram(a_data_yeast,0.2:0.01:0.7, ...
    'FaceColor',[233 113 50]/256,'EdgeColor',[1 1 1]);
xlim([0.2 0.7]);title("Distribution of a, S. cerevisiae",'FontWeight','Normal');
xlabel("a");ylabel("#Proteins");
ax = gca;
ax.LineWidth = 1;

%% Plot optimal a against beta and estimate beta from the median of yeast data
d_noLLPS = 0;
beta_vec = 0:0.01:10;
a_opt_human = beta_vec;
a_opt_yeast = beta_vec;
options = optimoptions('fmincon','display','none');
for i = 1:length(beta_vec)
    beta = beta_vec(i);
    a_opt_human(i) = fmincon(@(a)-payoff_overall(a,d_noLLPS,w_human,beta), ...
        0.5,[],[],[],[],0,1,[],options);
    a_opt_yeast(i) = fmincon(@(a)-payoff_overall(a,d_noLLPS,w_yeast,beta), ...
        0.5,[],[],[],[],0,1,[],options);
end
a_yeast_median = median(a_data_yeast);
beta_estimate = inv_linear_interpolation(beta_vec,a_opt_yeast,a_yeast_median);

figure;
subplot(1,2,1);
plot(beta_vec,a_opt_human,'LineWidth',3);
xlabel("\beta");ylabel("Optimal a under d=0");
title("H. sapiens",'FontWeight','Normal');
ax = gca;
ax.LineWidth = 1;
subplot(1,2,2);
plot(beta_vec,a_opt_yeast,'LineWidth',3);
xlabel("\beta");ylabel("Optimal a under d=0");
title("S. cerevisiae",'FontWeight','Normal');
hold on; scatter(beta_estimate,median(a_data_yeast),'Marker','hexagram');
ax = gca;
ax.LineWidth = 1;

%% Plot payoff function against parameter values
% Define range of parameter d and variable a
d = 0:1:3;
a = 0.001:0.001:0.999;
% Define the colormap
colors_spectrum = [0.84	0.10	0.11;
                    0.99	0.68	0.38;
                    0.67	0.87	0.64;
                    0.17	0.51	0.73];
figure;
subplot(2,2,1);
hold on;
for i = 1:length(d)
    f = payoff_overall(a,d(i),w_human,beta_estimate);
    plot(a,f,'Color',colors_spectrum(i,:),'LineWidth',3);
end
box on;
legend(["d=0","d=1","d=2","d=3"]);legend boxoff;
title("P_{max}(a), H. sapiens",'FontWeight','Normal');
xlabel("a");
ylabel("P_{max}(a)");
ax = gca;
ax.LineWidth = 1;

subplot(2,2,2);
hold on;
for i = 1:length(d)
    f = payoff_overall(a,d(i),w_yeast,beta_estimate);
    plot(a,f,'Color',colors_spectrum(i,:),'LineWidth',3);
end
box on;
legend(["d=0","d=1","d=2","d=3"]);legend boxoff;
title("P_{max}(a), S. cerevisiae",'FontWeight','Normal');
xlabel("a");
ylabel("P_{max}(a)");
ax = gca;
ax.LineWidth = 1;

subplot(2,2,3);
h = shannon_entropy20(a);
plot(a,h,'LineWidth',3);
title("Information capacity",'FontWeight','Normal');
xlabel("a");
ylabel("P_{info}(a)");
ax = gca;
ax.LineWidth = 1;

subplot(2,2,4);
A_human = optimal_abundance(a,w_human);
A_yeast = optimal_abundance(a,w_yeast);
plot(a,A_human,'LineWidth',3);
hold on;
plot(a,A_yeast,'LineWidth',3);
legend(["H. sapiens","S. cerevisiae"]);legend boxoff;
title("Maximal protein abundance",'FontWeight','Normal');
xlabel("a");
ylabel("A_{max}(a)");
ax = gca;
ax.LineWidth = 1;

%% Plot optimal a against d for human and yeast
d_vec = 0:0.01:30;
a_opt_human_varying_d = d_vec;
a_opt_yeast_varying_d = d_vec;
options = optimoptions('fmincon','display','none');
for i = 1:length(d_vec)
    a_opt_human_varying_d(i) = fmincon(@(a)-payoff_overall(a,d_vec(i),w_human, ...
        beta_estimate),0.5,[],[],[],[],0,1,[],options);
    a_opt_yeast_varying_d(i) = fmincon(@(a)-payoff_overall(a,d_vec(i),w_yeast, ...
        beta_estimate),0.5,[],[],[],[],0,1,[],options);
end
figure;
plot(d_vec,a_opt_human_varying_d,'LineWidth',3);
hold on;
plot(d_vec,a_opt_yeast_varying_d,'LineWidth',3);
xlabel("d");ylabel("a_{opt}");
title("Optimal a under varying d, human",'FontWeight','Normal');
ax = gca;
ax.LineWidth = 1;
legend(["H. sapiens","S. cerevisiae"]);legend boxoff;

%% Estimate d for each human protein
d_human = a_data_human;
for i = 1:length(d_human)
    d_human(i) = inv_linear_interpolation(d_vec, a_opt_human_varying_d, a_data_human(i));
end
figure;
histogram(d_human,50,'FaceColor',[233 113 50]/256,'EdgeColor',[1 1 1]);
xlabel("d");
ylabel("#Proteins");
title("Infered d for human proteins",'FontWeight','Normal');
ax = gca;
ax.LineWidth = 1;

%% Compute a under yeast w and human d
d_human_feasible = d_human(~isnan(d_human));
a_yeast_w_human_d = d_human_feasible;
for i = 1:length(a_yeast_w_human_d)
    a_yeast_w_human_d(i) = fmincon(@(a)-payoff_overall(a,d_human_feasible(i), ...
        w_yeast, beta_estimate),0.5,[],[],[],[],0,1,[],options);
end

figure;
histogram(a_yeast_w_human_d,50);
title("Optimal a under yeast w and human d",'FontWeight','Normal');
xlabel("a");
ylabel("#Proteins");
ax = gca;
ax.LineWidth = 1;

figure;
scatter(a_data_human(~isnan(d_human)),a_yeast_w_human_d,'Marker','*');
title("Optimal a under w of H. sapiens vs S. cerevisiae",'FontWeight','Normal');
xlabel("a, H.sapiens");
ylabel("a_{opt}, w of S. cerevisiae");
ax = gca;
ax.LineWidth = 1;
box on;

%% Definition of built-in functions
function f = payoff_func(a,d,w)
%Functional payoff function
f = (1 + d./(1 + exp(6.3 - 10.6*a)))./(1 + (w-1)*a);
end

function f = payoff_overall(a,d,w,beta)
%Overall payoff function
f = payoff_func(a,d,w) + beta*shannon_entropy20(a);
end

function h = shannon_entropy20(a)
%Information capacity term
h = -a.*log2(a/7) - (1-a).*log2((1-a)/13);
end

function A = optimal_abundance(a,w)
%Optimal abundance A under given a and w
A = 1./(w*a + 1-a);
end

function x0 = inv_linear_interpolation(xdata,ydata,y0)
% For a set of points (xdata,ydata) forming a curve and a known value y0,
% estimate a value x0 corresponding to y0 that falls on the curve
% The curve needs to be monotonic
delta_y = ydata - y0;
if min(delta_y) > 0
    x0 = NaN;
elseif max(delta_y) < 0
    x0 = NaN;
else
    deltay_sign_change_indicator = delta_y(1:end-1).*delta_y(2:end);
    idx_deltay_sign_change = find(deltay_sign_change_indicator < 0);
    xl = xdata(idx_deltay_sign_change);
    xr = xdata(idx_deltay_sign_change + 1);
    yl = ydata(idx_deltay_sign_change);
    yr = ydata(idx_deltay_sign_change + 1);
    x0 = xl + (xr - xl)*(y0 - yl)/(yr - yl);
end
end
    