clear all
clc

%% Script parameters
global ls h b Gcint Gcbulk Gceff E nu alpha1 alpha2
E = 126e9; % Young's modulus (Pa)
nu = 0.3; % Poisson's ratio
ls = 60e-6; % length scale (m)
Gcratio = 10; %Gcbulk/Gcint
b = 0.05e-3; % half interface thickness (m)
h = b/12; % interface element size (m)
Gcint = 562; % interface fracture toughness (N/m)
Gcbulk = Gcratio*Gcint; % bulk fracture toughness (N/m)
Gctol = 0.01; % Iterate until Gceff converges to within this tolerance
warning('on', 'all');

%% Yoshioka et al. (2021) as initial guess for ODE solver
Gceff = ((Gcbulk^2*h^2*exp((2*h)/ls) + Gcint^2*h^2*exp((2*h)/ls) +  Gcbulk^2*h^2*exp((4*b)/ls) + Gcint^2*h^2*exp((4*b)/ls) + 4*Gcbulk^2*ls^2*exp((2*h)/ls) + 4*Gcint^2*ls^2*exp((2*h)/ls) +  4*Gcbulk^2*ls^2*exp((4*b)/ls) + 4*Gcint^2*ls^2*exp((4*b)/ls) - 2*Gcbulk^2*h^2*exp((h + 2*b)/ls) + 2*Gcint^2*h^2*exp((h + 2*b)/ls) + 8*Gcbulk^2*ls^2*exp((h + 2*b)/ls) + 8*Gcint^2*ls^2*exp((h + 2*b)/ls) - 2*Gcbulk*Gcint*h^2*exp((2*h)/ls) + 2*Gcbulk*Gcint*h^2*exp((4*b)/ls) + 8*Gcbulk*Gcint*ls^2*exp((2*h)/ls) + 8*Gcbulk*Gcint*ls^2*exp((4*b)/ls) - 4*Gcbulk^2*h*ls*exp((2*h)/ls) + 4*Gcint^2*h*ls*exp((2*h)/ls) + 4*Gcbulk^2*h*ls*exp((4*b)/ls) + 4*Gcint^2*h*ls*exp((4*b)/ls) - 48*Gcbulk*Gcint*ls^2*exp((h + 2*b)/ls) + 8*Gcint^2*h*ls*exp((h + 2*b)/ls) +  8*Gcbulk*Gcint*h*ls*exp((4*b)/ls) - 24*Gcbulk*Gcint*h*ls*exp((h + 2*b)/ls))^(1/2) +  Gcbulk*h*exp(h/ls) + Gcint*h*exp(h/ls) - Gcbulk*h*exp((2*b)/ls) +  Gcint*h*exp((2*b)/ls) - 2*Gcbulk*ls*exp(h/ls) + 2*Gcint*ls*exp(h/ls) - 2*Gcbulk*ls*exp((2*b)/ls) + 2*Gcint*ls*exp((2*b)/ls)) /(2*h*exp(h/ls) + 2*h*exp((2*b)/ls) - 4*ls*exp(h/ls) + 4*ls*exp((2*b)/ls));
alpha1 = (exp((b - h/2)/ls)*(Gcbulk + Gceff)) / (Gcbulk*exp((b - h/2)/ls) - Gcbulk*exp(-(b - h/2)/ls) + Gceff*exp((b - h/2)/ls) + Gceff*exp(-(b - h/2)/ls));
alpha2 = (2*Gceff*exp((b - h/2)/ls)) / (Gcbulk*exp((b - h/2)/ls) - Gcbulk*exp(-(b - h/2)/ls) + Gceff*exp((b - h/2)/ls) + Gceff*exp(-(b - h/2)/ls));
fprintf('Yoshioka et al. (2021): Gceff =\t%.2f\n', Gceff);

%% Preliminary methods
fprintf('Solving for Gceff using preliminary methods...\n');
%total energy formulation (TEF) along orthogonal slice (Section 4.1)
[x_TEF, s_TEF] = iterativelySolve(@ftotal_orthog, @ftotal_orthog, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('total strain energy density formulation: Gceff =\t%.2f\n', Gceff)
%anisotropic strain energy density split (AES) along orthogonal slice (Section 4.2)
[x_AES, s_AES] = iterativelySolve(@fdev_orthog, @ftotal_orthog, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('Anisotropic strain energy density split: Gceff =\t%.2f\n', Gceff)

%% Maximum energy envelope method (Section 4.3)
fprintf('Generating maximum energy envelope SED functions...\n');
[kdev_handle, ktotal_handle] = SED_theta_max(nu, E, Gcint, ls);
global K_DEV_FUN K_TOTAL_FUN
K_DEV_FUN = kdev_handle;
K_TOTAL_FUN = ktotal_handle;
fprintf('Solving for Gceff using the maximum energy envelope method...\n');
[x, s] = iterativelySolve(@fdev, @ftotal, @bc, @guess_left, @guess_right, b, h, ls, Gcint, Gcbulk, Gctol);
fprintf('maximum energy envelope method: Gceff =\t%.2f\n', Gceff)

%% Data Extraction
s_data = [x, s];
H_raw = zeros(size(x));
idx_right = x >= 0;
k_vals_right = arrayfun(@(xi) K_TOTAL_FUN(Gcbulk, abs(xi)), x(idx_right)); 
H_raw(idx_right) = (ls * Gcbulk .* k_vals_right) ./ 2 ./ s(idx_right);
idx_left = ~idx_right;
k_vals_left = arrayfun(@(xi) K_DEV_FUN(Gcbulk, abs(xi)), x(idx_left));
H_raw(idx_left) = (ls * Gcbulk .* k_vals_left) ./ 2 ./ s(idx_left);
k_limit_val = K_TOTAL_FUN(Gceff, h/2); 
H_max_val = (ls * Gceff * k_limit_val) / 2;
H = min(H_raw, H_max_val);
H_data = [x, H];

%% Functions
function [x, s] = iterativelySolve(fdev, ftotal, bc, guess_left, guess_right, b, h, ls, Gcint, Gcbulk, Gctol)
    % The function f depends on Gceff, which is also being solved for. This
    % function therefore solves f iteratively, updating Gceff each time,
    % until Gceff is within tolerance.

    global Gceff

    xmesh_left = [linspace(0, b, 20000), linspace(b, 40*b, 20000)]; % starting from 0
    solinit_left = bvpinit(xmesh_left, guess_left);
    
    xmesh_right = [linspace(h/2, b, 20000), linspace(b, 40*b, 20000)];
    solinit_right = bvpinit(xmesh_right, guess_right);

    options = bvpset('RelTol', 1e-9, 'AbsTol', 1e-9);
    while true % iterate until Gceff is within tolerance
        sol_left = bvp4c(fdev, bc, solinit_left, options);
        sol_right = bvp4c(ftotal, bc, solinit_right, options);
        x = [-flipud(sol_left.x'); sol_right.x'];
        s = [flipud(sol_left.y(1,:)'); sol_right.y(1,:)'];
        DeltaGceff = Gceff_from_sProfile(x, s, b, ls, Gcint, Gcbulk) - Gceff; % change in Gceff
        Gceff = Gceff + DeltaGceff; % update Gceff
        if abs(DeltaGceff) < Gctol % is Gceff within tolerance?
            break;
        end
    end
end

function dsdx = ftotal_orthog(x, s, region)
    % The ODE to be solved from method 3.

    global ls Gceff Gcint Gcbulk nu h
    %Calculate k from sharp crack theory
    switch region
        case 1
            k = -(Gcint*(nu - 3))/(4*Gceff*ls*pi*abs(x)*s(1));
        case 2
            k = -(Gcint*(nu - 3))/(4*Gcbulk*ls*pi*abs(x)*s(1));
    end

    % Limit k to kmax
    kmax = -(Gcint*(nu - 3))/(4*Gceff*ls*pi*abs(h/2));
    if k > kmax
        k = kmax;
    end

    dsdx(1, 1) = s(2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function dsdx = fdev_orthog(x, s, region)
    % The ODE to be solved from method 3.

    global ls Gceff Gcint Gcbulk nu h
    %Calculate k from sharp crack theory
    switch region
        case 1
            k = (Gcint*(nu + 1))/(4*Gceff*ls*pi*abs(x)*s(1));
        case 2
            k = (Gcint*(nu + 1))/(4*Gcbulk*ls*pi*abs(x)*s(1));
    end

    % Limit k to kmax
    kmax = (Gcint*(nu + 1))/(4*Gceff*ls*pi*abs(h/2));
    if k > kmax
        k = kmax;
    end

    dsdx(1, 1) = s(2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function dsdx = fdev(x, s, region)
    global ls Gceff Gcbulk h K_DEV_FUN

    % Determine which Gc to use
    switch region
        case 1
            current_GcX = Gceff;
        case 2
            current_GcX = Gcbulk;
    end

    k = K_DEV_FUN(current_GcX, abs(x)) / s(1);
    kmax = K_DEV_FUN(Gceff, abs(h/2)); 
    
    % Limit k
    if k > kmax
        k = kmax;
    end
    
    dsdx(1, 1) = s(2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function dsdx = ftotal(x, s, region)
    global ls Gceff Gcbulk h K_TOTAL_FUN
    
    % Determine which Gc to use
    switch region
        case 1
            current_GcX = Gceff;
        case 2
            current_GcX = Gcbulk;
    end
    
    k = K_TOTAL_FUN(current_GcX, abs(x)) / s(1);
    kmax = K_TOTAL_FUN(Gceff, abs(h/2)); 
    
    % Limit k
    if k > kmax
        k = kmax;
    end
    
    dsdx(1, 1) = s(2);
    dsdx(2, 1) = s(1)*(k + 1/ls^2) - 1/ls^2;
end

function res = bc(YL, YR)
    % YL(i, j) is value of s or s' at left boundary of region j.
    % YR(i, j) is value of s or s' at right boundary of region j.
    % i = 1 for s, or i = 2 for s'.
    % For example, YR(2, 1) is s' on right of region 1.
    
    global Gcbulk Gceff
    res = [YL(1, 1); % s(0) = 0
        YR(1, 1) - YL(1, 2); % s1(b) = s2(b)
        Gceff*YR(2, 1) - Gcbulk*YL(2, 2); % first Weierstrass-Erdmann corner condition
        YR(1, 2) - 1]; % s(Inf) = 1;
end

function g = guess_left(x, ~)
    % Use Yoshioka et al.'s (2021) s profile without the element-size
    % effect as the initial guess for the left side.

    global alpha1 alpha2 ls
    if abs(x) > 0
        g = [1 - alpha2*exp(-(abs(x))/ls);
            (alpha2*exp((-abs(x))/ls))/ls];
    else
        g = [1 - alpha1*exp(-(abs(x))/ls) - (1 - alpha1)*exp((abs(x))/ls);
            (alpha1*exp((-abs(x))/ls))/ls + (exp(-(-abs(x))/ls)*(alpha1 - 1))/ls];
    end
end

function g = guess_right(x, ~)
    % Use Yoshioka et al.'s (2021) s profile with the element-size effect
    % as the initial guess for the right side.

    global alpha1 alpha2 h ls
    if abs(x) < h/2
        g = [0; 0];
    elseif abs(x) > h/2
        g = [1 - alpha2*exp(-(abs(x) - h/2)/ls);
            (alpha2*exp((h/2 - abs(x))/ls))/ls];
    else
        g = [1 - alpha1*exp(-(abs(x) - h/2)/ls) - (1 - alpha1)*exp((abs(x) - h/2)/ls);
            (alpha1*exp((h/2 - abs(x))/ls))/ls + (exp(-(h/2 - abs(x))/ls)*(alpha1 - 1))/ls];
    end
end

function Gceff = Gceff_from_sProfile(x, s, b, ls, Gcint, Gcbulk)
    % Apply the energy balance equation to compute Gceff
    B = unique([x s],'rows');
    x = B(:,1); s = B(:,2);

    bulk_left  = x < -b;
    intf  = abs(x) <= b;
    bulk_right = x > b;
    
    ds_dx = gradient(s)./gradient(x);
    S = 1/2*((1 - s).^2/ls + ls*ds_dx.^2);
    Gceff = (Gcint - Gcbulk*trapz(x(bulk_left), S(bulk_left)) - Gcbulk*trapz(x(bulk_right), S(bulk_right)))/trapz(x(intf), S(intf));
end

function [kdev_fun, ktotal_fun] = SED_theta_max(nu_val, E_val, Gcint_val, ls_val)
    %Compute the maximum SED functions for use in the governing equation.
    %k=2H(x)/ls/Gc(x)
    syms theta r GcX positive real
    KII_val = sqrt(E_val * Gcint_val);
    fxx_II = -sin(theta/2)*(2 + cos(theta/2)*cos(3*theta/2));
    fyy_II = sin(theta/2)*cos(theta/2)*cos(3*theta/2);
    fxy_II = cos(theta/2)*(1 - sin(theta/2)*sin(3*theta/2));
    
    sigxx = (1/sqrt(2*pi*r))*KII_val*fxx_II;
    sigyy = (1/sqrt(2*pi*r))*KII_val*fyy_II;
    sigxy = (1/sqrt(2*pi*r))*KII_val*fxy_II;

    epsxx = (1/E_val)*(sigxx - nu_val*sigyy);
    epsyy = (1/E_val)*(sigyy - nu_val*sigxx);
    epsxy = (1+nu_val)/E_val*sigxy;

    U = simplify(0.5*sigxx*epsxx + 0.5*sigyy*epsyy + sigxy*epsxy);

    sigvol = (sigxx + sigyy)/2;
    epsvol = sigvol/(E_val/2/(1 - nu_val));
    Uvol = 0.5*sigvol*epsvol;
    Udev = simplify(U - Uvol);

    % Optimisation
    % substitute r = 1/sin(theta) for finding the angle
    U_scan = subs(U, r, 1/sin(theta));
    Udev_scan = subs(Udev, r, 1/sin(theta));
    
    f_opt = matlabFunction(-U_scan, 'Vars', theta);
    fdev_opt = matlabFunction(-Udev_scan, 'Vars', theta);
    
    theta_crit = fminbnd(f_opt, 0, pi);
    thetadev_crit = fminbnd(fdev_opt, 0, pi);
    
    syms x_dist positive real % Create a clean distance variable
    
    % Substitute the geometric relation: r = dist_from_centreline/sin(theta)
    U_final = subs(U, {theta, r}, {theta_crit, x_dist/sin(theta_crit)});
    ktotal_sym = simplify(2 * U_final / (ls_val * GcX));
    
    Udev_final = subs(Udev, {theta, r}, {thetadev_crit, x_dist/sin(thetadev_crit)});
    kdev_sym = simplify(2 * Udev_final / (ls_val * GcX));

    % Generate Handles: k = f(GcX, dist)
    ktotal_fun = matlabFunction(ktotal_sym, 'Vars', {GcX, x_dist});
    kdev_fun   = matlabFunction(kdev_sym,   'Vars', {GcX, x_dist});
end
