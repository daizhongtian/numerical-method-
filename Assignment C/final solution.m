


fprintf('\n=== start Task 1: demonstrate ODE  ===\n');

tDemo = 0:0.01:5;
x0_demo = 450;  y0_demo = 50;
alpha_demo = [10, -0.1, -0.004];
beta_demo  = [-10, 0.06, -0.001];
[x1, y1] = predatorPreySolve('ode45', tDemo, x0_demo, y0_demo, alpha_demo, beta_demo);
[x2, y2] = predatorPreySolve('gear2', tDemo, x0_demo, y0_demo, alpha_demo, beta_demo);
[x3, y3] = predatorPreySolve('kutta', tDemo, x0_demo, y0_demo, alpha_demo, beta_demo);
figure('Name','Task 1 – ODE45  (blue: x, red: y)');
plot(tDemo, x1, 'b-', tDemo, y1, 'r-','LineWidth',1.5);
xlabel('t'); ylabel('Population'); title('Task 1 – ode45 ');

figure('Name','Task 1 – Gear2  (blue: x, red: y)');
plot(tDemo, x2, 'b--', tDemo, y2, 'r--','LineWidth',1.5);
xlabel('t'); ylabel('Population'); title('Task 1 – Gear2 ');

figure('Name','Task 1 – RK4(Kutta)  (blue: x, red: y)');
plot(tDemo, x3, 'b:', tDemo, y3, 'r:','LineWidth',1.5);
xlabel('t'); ylabel('Population'); title('Task 1 – RK4 (Kutta) ');

fprintf('\n=== above is Task 1 demonstration, now proceeding to Task 2 ===\n');

[aH, bH, info2] = predatorPreyEstimate('animals_1.csv');
fprintf('\nTask 2 completed: parameter estimation results are as follows:\n');
fprintf('  alphaHat = [%g  %g  %g]\n', aH);
fprintf('  betaHat  = [%g  %g  %g]\n', bH);

fprintf('\n=== proceeding to Task 3 ===\n');

results = predatorPreyTask3('Data','animals_1.csv', ...
            'Tol' , [1e-2 1e-3 1e-4 1e-6 1e-8], ...
            'Step', [0.04 0.02 0.01 0.005 0.0025]);

fprintf('\nAll tasks have been completed.\n');



function [xHat, yHat] = predatorPreySolve(method, tVec, x0, y0, alpha, beta, varargin)


    if nargin == 0
        fprintf('=== no parameter ===\n');
        predatorPreyEstimate('animals_1.csv');
        return;
    end



    tVec = tVec(:);
    switch lower(string(method))
        case "ode45"
            [xHat, yHat] = solverOde45(tVec, x0, y0, alpha, beta, varargin{:});
        case "gear2"
            [xHat, yHat] = solverGear2(tVec, x0, y0, alpha, beta);
        case "kutta"
            [xHat, yHat] = solverKutta(tVec, x0, y0, alpha, beta);
        otherwise
            error('Unknown method "%s". select "ode45","gear2" or "kutta".', method);
    end
end 





function [xHat, yHat] = solverOde45(t, x0, y0, alpha, beta, varargin)


    u0  = [x0; y0];
    rhs = @(~, u) rhsPredPrey(u, alpha, beta);

    if isempty(varargin)
        sol = ode45(rhs, [t(1), t(end)], u0);
    else
        sol = ode45(rhs, [t(1), t(end)], u0, odeset(varargin{:}));
    end

    U    = deval(sol, t).';  
    xHat = U(:, 1);
    yHat = U(:, 2);
end



function [xHat, yHat] = solverGear2(t, x0, y0, alpha, beta)


    N = numel(t);
    if N < 2
        error('Gear2 requires at least two time points.');
    end

    h_all = diff(t);
    tol_h = eps(max(t)) * 100;
    if any(abs(h_all - h_all(1)) > tol_h)
        error('solverGear2 only supports a uniform grid; detected non-uniform t.');
    end
    h = h_all(1);

    xHat = zeros(N, 1);
    yHat = zeros(N, 1);
    xHat(1) = x0;
    yHat(1) = y0;

    u      = [x0; y0];
    f_prev = rhsPredPrey(u, alpha, beta); 

    % Heun (RK2) 
    u_tilde = u + h * f_prev;
    f_tilde = rhsPredPrey(u_tilde, alpha, beta);
    u       = u + (h/2) * (f_prev + f_tilde);
    xHat(2) = u(1);
    yHat(2) = u(2);

    %  Adams–Bashforth 2 
    for n = 2 : N-1
        f_curr = rhsPredPrey(u, alpha, beta);
        u      = u + h * (1.5 * f_curr - 0.5 * f_prev);
        f_prev = f_curr;
        xHat(n+1) = u(1);
        yHat(n+1) = u(2);
    end
end


%----------------------------------------------------------------------

function [xHat, yHat] = solverKutta(t, x0, y0, alpha, beta)


    N    = numel(t);
    xHat = zeros(N, 1);
    yHat = zeros(N, 1);
    xHat(1) = x0;
    yHat(1) = y0;

    u = [x0; y0];
    for n = 1 : (N-1)
        h = t(n+1) - t(n);

        k1 = rhsPredPrey(u,              alpha, beta);
        k2 = rhsPredPrey(u + (h/2)*k1,   alpha, beta);
        k3 = rhsPredPrey(u + (h/2)*k2,   alpha, beta);
        k4 = rhsPredPrey(u +     h *k3,  alpha, beta);

        u = u + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        xHat(n+1) = u(1);
        yHat(n+1) = u(2);
    end
end


%----------------------------------------------------------------------

function dudt = rhsPredPrey(u, alpha, beta)


    x = u(1);
    y = u(2);

    
    dudt = [ alpha(1)*x + alpha(2)*x.*y + alpha(3)*x.^2;
             beta(1)*y  + beta(2)*x.*y  + beta(3)*y.^2 ];
end
function [alphaHat, betaHat, info] = predatorPreyEstimate( ...
                           dataFile, initGuess, odeOpts, fminOpts)

if nargin<2 || isempty(initGuess)
    initGuess = [10 -0.03 -0.001  -8 0.02 -0.004];
end
if nargin<3 || isempty(odeOpts)
    odeOpts = {'RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-2};
end
if nargin<4 || isempty(fminOpts)
    fminOpts = optimset('Display','iter', ...
                        'TolX',1e-6,'TolFun',1e-6, ...
                        'MaxIter',4e3,'MaxFunEvals',4e3);
end


M = readmatrix(dataFile);
if size(M,2) < 3
    error('Data file must have at least 3 columns: [t,xObs,yObs].');
end
tAll = M(:,1);
xObs = M(:,2);
yObs = M(:,3);


x0 = xObs(1);
y0 = yObs(1);


theta0 = [ log(abs(initGuess(1)));          % θ1 = log(a1)
           initGuess(2:3).';               % θ2,θ3 = a2,a3
           log(abs(initGuess(4)));          % θ4 = log(|b1|)
           initGuess(5:6).' ];             % θ5,θ6 = b2,b3

theta2param = @(th)[ ...
     exp(th(1)),  th(2),  th(3), ...   % a1 a2 a3
    -exp(th(4)),  th(5),  th(6) ];     % b1 b2 b3   (b1<0)


    function rss = objFun(theta)       
        p     = theta2param(theta);
        alpha = p(1:3);   beta  = p(4:6);

        % integrate with ode45
        try
            rhs = @(t,u)[ alpha(1)*u(1)+alpha(2)*u(1)*u(2)+alpha(3)*u(1).^2 ; ...
                           beta(1)*u(2)+beta(2)*u(1)*u(2)+beta(3)*u(2).^2 ];
            sol = ode45(rhs,[tAll(1) tAll(end)],[x0; y0],odeset(odeOpts{:}));
            xyFit = deval(sol,tAll).';        % N×2
        catch
            rss = 1e20;   
            return
        end

        dx = xyFit(:,1) - xObs;
        dy = xyFit(:,2) - yObs;
        rss = sum(dx.^2 + dy.^2);            
    end

[thetaHat,fval,exitflag,output] = fminsearch(@objFun,theta0,fminOpts);

pHat      = theta2param(thetaHat);
alphaHat  = pHat(1:3).';        
betaHat   = pHat(4:6).';
rss       = fval;

fprintf('\n=== Task 2 finished (solver=ode45, search=fminsearch) ===\n');
fprintf('alphaHat = [%g  %g  %g]\n', alphaHat);
fprintf('betaHat  = [%g  %g  %g]\n', betaHat);
fprintf('RSS      = %.6g\n', rss);

try
    rhs = @(t,u)[ alphaHat(1)*u(1)+alphaHat(2)*u(1)*u(2)+alphaHat(3)*u(1).^2 ; ...
                  betaHat(1)*u(2)+betaHat(2)*u(1)*u(2)+betaHat(3)*u(2).^2 ];
    sol = ode45(rhs,[tAll(1) tAll(end)],[x0; y0],odeset(odeOpts{:}));
    xyFit = deval(sol,tAll).';
    figure('Name','Task 2 – Observation vs Model');
    plot(tAll,xObs,'bo',tAll,yObs,'ro','MarkerFaceColor','auto'); hold on
    plot(tAll,xyFit(:,1),'b-','LineWidth',1.4);
    plot(tAll,xyFit(:,2),'r-','LineWidth',1.4);
    legend('x obs','y obs','x model','y model','Location','best');
    xlabel('t'); ylabel('Population'); grid on
    title('Parameter estimation (ode45 + fminsearch)');
catch ME
    warning(ME.identifier,'Plotting failed: %s',ME.message);
end

info = struct('rss',rss,'exitflag',exitflag,'output',output, ...
              'initGuess',initGuess,'odeOpts',{odeOpts});
end



% Task-3  


function results = predatorPreyTask3(varargin)

p = inputParser;
addParameter(p,'Data', 'animals_1.csv', @ischar);

addParameter(p,'Tol',  [1e-2 1e-3 1e-4 1e-6 1e-8], ...
             @(v)isnumeric(v)&&isvector(v));

addParameter(p,'Step', [ ...
    0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04, ...
    0.06, 0.08, 0.1, 0.12, 0.16, 0.2, ...
    0.3, 0.4, 0.6, 0.8, 1.0 ], @(v)isnumeric(v)&&isvector(v));


parse(p,varargin{:});
dataFile = p.Results.Data;
tolList  = p.Results.Tol(:).';
hList    = p.Results.Step(:).';

M    = readmatrix(dataFile);
tObs = M(:,1);  xObs = M(:,2);  yObs = M(:,3);
x0 = xObs(1);   y0 = yObs(1);
dtMin = min(diff(tObs));

pTrue = getTrueParams(dataFile);   

cfg = struct('label',{},'method',{},'tol',{},'h',{},...
             'rmsA',{},'rmsB',{},'wallTime',{},'RE',{},'exitFlag',{});

for tol = tolList
    cfg(end+1) = struct('label',sprintf('ode45 Tol=%g',tol), ...
        'method',"ode45",'tol',tol,'h',NaN,'rmsA',NaN,'rmsB',NaN, ...
        'wallTime',NaN,'RE',[],'exitFlag',0);
end
for h = hList
    cfg(end+1) = struct('label',sprintf('gear2 h=%g',h), ...
        'method',"gear2",'tol',NaN,'h',h,'rmsA',NaN,'rmsB',NaN, ...
        'wallTime',NaN,'RE',[],'exitFlag',0);
end
for h = hList
    cfg(end+1) = struct('label',sprintf('rk4   h=%g',h), ...
        'method',"kutta",'tol',NaN,'h',h,'rmsA',NaN,'rmsB',NaN, ...
        'wallTime',NaN,'RE',[],'exitFlag',0);
end

fprintf('\n===========  Task-3   ===========\n');
Nrep = 10;                        
for k = 1:numel(cfg)
    fprintf('>> %-18s : ', cfg(k).label);
    
   
    times = zeros(Nrep,1);
    for j = 1:Nrep
        t0 = tic;
        localFit(cfg(k), tObs, xObs, yObs, x0, y0, dtMin, pTrue);
        times(j) = toc(t0);
    end
    cfg(k).wallTime = mean(times);
    
   
    [rmsA, rmsB, reVec, flag] = ...
        localFit(cfg(k), tObs, xObs, yObs, x0, y0, dtMin, pTrue);
    cfg(k).rmsA     = rmsA;
    cfg(k).rmsB     = rmsB;
    cfg(k).RE       = reVec;
    cfg(k).exitFlag = flag;
    
    fprintf('avg time %.4gs   RMS-α %.3g   RMS-β %.3g\n', ...
            cfg(k).wallTime, cfg(k).rmsA, cfg(k).rmsB);
end
results = cfg;
%—— 4. draw ——----------------------------------------------------------
plotTask3Results(results);
end   




function [rmsA,rmsB,reVec,flag] = localFit(c,tObs,xObs,yObs,...
                                           x0,y0,dtMin,pTrue)
flag = 0;


theta0 = [log(10); -0.03; -1e-3; log(8); 0.02; -4e-3];
mapP   = @(th)[exp(th(1)), th(2), th(3), -exp(th(4)), th(5), th(6)];

    function rss = obj(th)
        p = mapP(th);  a=p(1:3);  b=p(4:6);
        try
            if c.method=="ode45"
                [xs,ys] = predatorPreySolve('ode45', tObs,...
                    x0,y0,a,b, 'AbsTol',c.tol,'RelTol',c.tol);
            else
    tInt = tObs(1):c.h:tObs(end);
    if tInt(end) < tObs(end), tInt = [tInt, tObs(end)]; end
    [xt, yt] = predatorPreySolve(c.method, tInt, x0, y0, a, b);

   
    if any(isnan(xt)) || any(isnan(yt))
        rss = 1e20;
        return;
    end

    xs = interp1(tInt, xt, tObs, 'pchip', 'extrap');
    ys = interp1(tInt, yt, tObs, 'pchip', 'extrap');
end
        catch
            rss = 1e20; return;
        end
        rss = sum((xs-xObs).^2 + (ys-yObs).^2);  
    end

opts = optimset('Display','off','MaxIter',1e3,'MaxFunEvals',1e3,...
                'TolX',1e-6,'TolFun',1e-6);
thetaHat = fminsearch(@obj,theta0,opts);


pEst = mapP(thetaHat);          
reVec= abs(pEst - pTrue) ./ abs(pTrue);   

rmsA = sqrt(mean(reVec(1:3).^2));
rmsB = sqrt(mean(reVec(4:6).^2));
end   



function plotTask3Results(cfg)
fprintf('\n>>> draw Task-3 curve …\n');
T = struct2table(cfg);

valid = T.rmsA < Inf & T.rmsB < Inf;
T = T(valid,:);

isO45 = startsWith(T.label,"ode45");
isG2  = startsWith(T.label,"gear2");
isRK4 = startsWith(T.label,"rk4");

%—— Fig-A  ode45 : Tol → precision  -----------------------------------------
figure('Name','ode45 – Tol vs Parameter Accuracy','Color','w');
tiledlayout(2,1,'TileSpacing','compact');

nexttile
semilogx(T.tol(isO45), T.rmsA(isO45),'bo-','LineWidth',1.6);
ylabel('RMS RE_{α}'); grid on;
title('ode45 :  Tol  → α-precision');

nexttile
semilogx(T.tol(isO45), T.rmsB(isO45),'rs-','LineWidth',1.6);
ylabel('RMS RE_{β}'); xlabel('Tol (RelTol = AbsTol)'); grid on;
title('ode45 :  Tol  → β-precision');

%—— Fig-B  step : h → precision  ----------------------------------------
figure('Name','Fixed-step – h vs Parameter Accuracy','Color','w');
tiledlayout(2,1,'TileSpacing','compact');

nexttile
loglog(T.h(isG2), T.rmsA(isG2),'r^-','LineWidth',1.6); hold on;
loglog(T.h(isRK4),T.rmsA(isRK4),'gs-','LineWidth',1.6);
ylabel('RMS RE_{α}'); grid on;
title('Fixed-step :  h  → α-precision');
legend('gear2','rk4','Location','best');

nexttile
loglog(T.h(isG2), T.rmsB(isG2),'r^-','LineWidth',1.6); hold on;
loglog(T.h(isRK4),T.rmsB(isRK4),'gs-','LineWidth',1.6);
ylabel('RMS RE_{β}'); xlabel('Step size  h'); grid on;
title('Fixed-step :  h  → β-precision');

%—— Fig-C :  ode45   Tol → wall-time ----------------------------------
figure('Name','ode45 – Tol vs Wall-clock Time','Color','w');
semilogx(T.tol(isO45), T.wallTime(isO45), 'bo-', 'LineWidth', 1.6);
ylabel('Wall-clock time [s]');
xlabel('Tol  (RelTol = AbsTol)');
title('ode45 :  Tol  →  speed');
grid on;

%—— Fig-D :  fixed-step  h → wall-time --------------------------------
figure('Name','Fixed-step – h vs Wall-clock Time','Color','w');
loglog(T.h(isG2 ),  T.wallTime(isG2 ),  'r^-', 'LineWidth', 1.6); hold on;
loglog(T.h(isRK4), T.wallTime(isRK4), 'gs-', 'LineWidth', 1.6);
ylabel('Wall-clock time [s]');
xlabel('Step size  h');
title('Fixed-step :  h  →  speed');
legend('gear2', 'rk4', 'Location', 'best');
grid on;

end   




function pTrue = getTrueParams(fname)
switch fname
    case 'animals_1.csv'
        pTrue = [10 -0.03 -0.001  -8 0.02 -0.004];
    case 'animals_2.csv'
        pTrue = [11 -0.04 -0.002 -13 0.03 -0.003];
    case 'animals_3.csv'
        pTrue = [7  -0.01 -0.003 -15 0.01 -0.005];
    otherwise
        error('Unknown data file: %s', fname);
end
end
