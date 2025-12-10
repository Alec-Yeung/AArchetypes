%% Grid over w and d, solve optimal a,A for each (w,d), and plot filled contours of a
% ----- ranges -----
w_vec = linspace(0, 6, 100);       
d_vec = linspace(0, 6, 100);       
nw = numel(w_vec);
nd = numel(d_vec);
c_ext = 1;

% ----- storage -----
a_grid   = zeros(nd, nw);         
A_grid   = zeros(nd, nw);
obj_grid = zeros(nd, nw);
ext_grid = zeros(nd, nw);

% ----- bounds & options -----
lb = [0; 0];
ub = [1; Inf];
opts = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','sqp-legacy', ...
    'OptimalityTolerance',1e-9, ...
    'StepTolerance',1e-9, ...
    'ConstraintTolerance',1e-9);

% ----- solve for each (w,d) -----
for iw = 1:nw
    w = w_vec(iw);
    for id = 1:nd
        d = d_vec(id);

        fun     = @(x) benefit_Human(x, d);
        nonlcon = @(x) budget_constraint_Human(x, w);

        x0 = [0.5; 0.5];

        [xbest, fval] = fmincon(fun, x0, [], [], [], [], lb, ub, nonlcon, opts);

        a_grid(id, iw)   = xbest(1);
        A_grid(id, iw)   = xbest(2);
        obj_grid(id, iw) = -fval;
        ext_grid(id, iw) = external_cost_Human(xbest(1), c_ext);
    end
end

% ----- plot1: filled contour of a over (w,d) -----
[W, D] = meshgrid(w_vec, d_vec);

figure; 
contourf(W, D, a_grid, 30, 'LineColor', 'none');
axis tight;
set(gca, 'YDir', 'normal'); 
grid off;
xlabel('w'); 
ylabel('d'); 
title('Optimal a across (w, d)');
c_low  = [158, 202, 254] / 255;   % #9ecafe
c_high = [254, 206, 230] / 255;   % #fecee6
N = 256;
cmap = [linspace(c_low(1),  c_high(1),  N)', ...
        linspace(c_low(2),  c_high(2),  N)', ...
        linspace(c_low(3),  c_high(3),  N)'];
colormap(cmap);
caxis([0 1]);                     
cb = colorbar; cb.Label.String = 'optimal a';
hold on;
xline(5.77, '--k', 'LineWidth', 1);  
xline(0.26, '--k', 'LineWidth', 1);  
xline(0.54, '--k', 'LineWidth', 1);  
hold off;

% ----- plot2: filled contour of external cost over (w,d) -----
figure;
contourf(W, D, a_grid, 30, 'LineColor', 'none');  
axis tight; 
set(gca, 'YDir', 'normal'); 
grid off;
xlabel('w'); 
ylabel('d'); 
title('External cost across (w, d)');
colormap(cmap);
caxis([0 c_ext]);
cb = colorbar; cb.Label.String = 'external cost (c \cdot p_A)';
hold on;
xline(5.77, '--k', 'LineWidth', 1);  
xline(0.26, '--k', 'LineWidth', 1);  
xline(0.54, '--k', 'LineWidth', 1);  
hold off;


%%
function p = benefit_Human(x,d)
a = x(1);
A = x(2);
pA = LLPS_prob_Human(a);
benefit = A*(1 + d*pA);
p = -benefit;
end

function p = benefit_Yeast(x,d)
a = x(1);
A = x(2);
pA = LLPS_prob_Yeast(a);
benefit = A*(1 + d*pA);
p = -benefit;
end

function p = LLPS_prob_Human(x)
p = 1 / (1 + exp(6.3296 - 10.5557*x));
end

function p = LLPS_prob_Yeast(x)
p = 1 / (1 + exp(7.0094 - 11.5778*x));
end

function [c,ceq] = budget_constraint_Human(x,w)
a = x(1);
A = x(2);
c = (w*a + (1-a)) * A - 1;
ceq = [];
end

function [c,ceq] = budget_constraint_Yeast(x,w)
a = x(1);
A = x(2);
c = (w*a + (1-a)) * A - 1;
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
