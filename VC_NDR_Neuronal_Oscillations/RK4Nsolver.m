function [sol,allstates,max_err] = RK4Nsolver(lims,h_vec,out1,prev,images,state,params)

% Function "RK4Nsolver.m" v3.0.0, tested 29 June 2021
% Written by T.D. Brown Fall 2019
%
% THIS FUNCTION DEPENDS ON DERIVFUNC.M!!
% 
% This function is a numerical solver and solves the N-DIMENSIONAL 
% non-linear system of differential equations, { x1'=f1(t,x1,x2,...,xn),
% x2'=f2(t,x1,x2,...,xn), ..., xn'=fn(t,x1,x2,...,xn) } using
% the Runge-Kutta 4th order numerical method, which is more robust to
% numerical instability and has order of convergence of 4 (the error 
% between subsequent approximations decreases as the 4th power of the time
% step). This function has some nice features for investigating convergence
% and allowing user-defined time-stepping. 
% It also plots the solutions, if desired.
%
% Recent updates have incorporated local state-dependence which is carried
% through between calls to derivfunc, useful for e.g., incorporating
% inherent hysteresis models
%
% The ability to work for any generalized RK method will be added in a
% later version. For now, the function uses the "standard 1/3" rule. 
%
% This matlab function will return an array of time t, x1(t), x2(t), ... 
% xn(t), and an error parameter, relative to a previous run of the solver.
% This is useful, for example, when studying convergence as f(timestepping)
%
% Output: sol = [time_vec; x1_vec; x2_vec; xn_vec], [array with state variables in ROWS, timesteps in COLUMNS]
% Output: err = log10([max_error_x / avg_val_x, max_error_y / avg_val_y]), [ROW vector]
% Input: lims = [initial_time,final_time], [ROW vector]
% Input: h_vec = [delta_t1, delta_t2, deltat_n-1], [ROW vector] OR
% Input: h_vec = Length of uniform step sizes between lims [scalar]
% Input: out1 = [x1(initial_time);x2(initial_time);...xn(iniial_time)], [COLUMN vector]
% Input: prev = sol array from previous run, [array]
% Input: imgflag = [desired_var1, ..., desired_vark], [ROW vector of natural numbers]
% Figure 1: x(desired_var1)(t) for current and previous solution
% Figure 2: x(desired_var2)(t) for current and previous solution
% ...
% Figure k: x(desired_vark)(t) for current and previous solutions
% Figure k+1: Difference plot between current and previous solution for k
% desired variables
%
% Example 1: Solve the system x'=y, y'=-by-kx from t=0 s to t=10 s, using
% 100 evenly spaced steps, subject to ICs, x(t0)=1, y(t0)=0
% >> lims = [0,10];
% >> h_vec = 0.01; % Gives 100 evenly spaced time-steps
% >> out1 = [1 0]' (Subject to ICs)
% >> prev = []; % We don't have a previous solution yet
% >> [sol100,~,err] = RK4n2solver(lims,h_vec,out1,prev);
% >> ***DERIVFUNC.M must be appropriately programmed for x'=y, y'=-by-kx***
%
% Example 2: Solve the same system, but now compare the solutions
% when we increase the timesteps from 100 to 1,000
% Notice how the solution takes longer but is better converged
% >> h_vec = 0.001;
% >> prev = sol100;
% >> [sol1000,~,err] = RK4n2solver(lims,h_vec,out1,prev);
%
% Example 3: Try increasing the number of timesteps to 10,000
% The relative error in x(t) and y(t) for 10,000 steps compared to 1,000 is
% about 10% (10^-1.02), so our solution with 1,000 steps has converged within only 10%
% Similarly the solution with 10,000 steps has converged within 1% (as you can check)
% >> h_vec = 0.0001;
% >> prev = sol1000
% >> [sol10000,err] = RK4n2solver(lims,h_vec,out1,prev)
%
% Example 4: Now try a non-uniform time-step spacing that is bunched
% together at the beginning but then spreads out
% >> h_vec = [5:-0.05:1]; % 80 non-evenly spaced points
% >> prev = [];
% >> [sol,err] = RK4n2solver(lims,h_vec,out1,prev);
%
% To.do: Implement functionality for general RK methods
% To.do: Play with ability to make plots (choose among N variables)

%% User-Modified Block / Internal Parameters
% There aren't any! Change function behavior through inputs only

%% Don't Change Anything Below This!
close all

%% Runge-Kutta Solver
% This uses lims and h_vec to define the vector of time points, "in"
if size(h_vec,2) == 1
    lgt = round(1/h_vec);
    h_vec = 1/lgt*ones(1,lgt);
end

if sum(h_vec) ~= 1
    h_vec = h_vec / sum(h_vec);
end

in = lims(1)+(lims(2)-lims(1))*[0,cumsum(h_vec)];

% Initializes output variables
out = zeros(size(out1,1),size(in,2));
[~,state0] = derivfunc(in(1),out(:,1),state,params);
allstates(size(in,2)) = struct();

for fn = fieldnames(state0)'
    allstates(1).(fn{1}) = state0.(fn{1});
end

if ~any(size(prev))
   prev  = [in;out];
end

out(:,1) = out1;

% Actually does the numerical solving, using RK4 method
% Notice the crucial role of derivfunc in defining the problem
coeff = [ [0, 0]', [1/2, 1/2]', [1/2 1/2]', [1, 1]' ]; 
kzero = zeros(size(out1,1),size(coeff,2)+1); 
prog = 0;
    for ii = 1:length(in)-1
        cur_h = in(ii+1)-in(ii);
        cur_t = in(ii);
        cur = out(:,ii);    
        k = kzero;  
        
        % Don't Update State Until End of RK4 Execution
        for jj = 1:size(k,2)-1
            [output,~] = derivfunc(cur_t+coeff(1,jj)*cur_h, ...
                cur+coeff(2,jj)*k(:,jj),state,params);
            k(:,jj+1) = cur_h*output;
        end
            
        nxt = cur + (k*[0 1/6 2/6 2/6 1/6]'); % This zero is very important!
        nxt_t = in(ii+1);
        out(:,ii+1) = nxt;
        
        % Only update the state after all function calls for single step of 
        % RK4 solver are completed!
        [~,state] = derivfunc(nxt_t,nxt,state,params);
        allstates(ii+1) = state;
        
        % Break if a NaN is detected
        if any(isnan(nxt))
            disp('Nan detected, failed to converge!')
            break
        end     
        
        % Read out continuous simulation progress by every 1%
        if 100*ii/length(h_vec)>prog
        	disp(['Progress: ',num2str(prog),'%']);
            prog = prog+1;
        end
        
    end

% Integration Complete, Return Solution Array    
sol = [in;out];
    
%% Interpolate to allow comparisons of current and previous runs
if size(sol,2)>=size(prev,2) % New solution is higher resolved than old
    in2 = in;  
    inn2 = prev(1,:);
    x2 = out;
    xx1 = prev(2:size(prev,1),:);
    x1 = zeros(size(out,1),length(in2));
    xx2 = zeros(size(out,1),length(inn2));
    for kk = 1:size(out,1)
        x1(kk,:) = interp1(prev(1,:),prev(kk+1,:),in2);
        xx2(kk,:)=interp1(in,out(kk,:),inn2);
    end
else                         % Old solution is lower resolved than old
    in2 = prev(1,:);
    inn2 = in;
    x1 = prev(2:size(prev,1),:);
    xx2 = out;
    x2 = zeros(size(out,1),length(in2));
    xx1 = zeros(size(out,1),length(inn2));
    for kk = 1:size(out,1)
        x2(kk,:) = interp1(in,out(kk,:),in2);
        xx1(kk,:) = interp1(prev(1,:),prev(kk+1,:),inn2);
    end
end

% Because MATLab interp1 function pads extrapolation with nan...
x1(:,1) = prev(2:size(prev,1),1); % First column
x1(:,size(x1,2)) = prev(2:size(prev,1),size(prev,2)); % Last column
x2(:,1) = out(:,1); % First column
x2(:,size(x2,2)) = out(:,size(out,2)); % Last column
xx1(:,1) = prev(2:size(prev,1),1); % First column
xx1(:,size(xx1,2)) = prev(2:size(prev,1),size(prev,2)); % Last column
xx2(:,1) = out(:,1); % First column
xx2(:,size(xx2,2)) = out(:,size(out,2)); % Last column

%% Calculating Error Parameters
% Error parameter: log10 of maximum of abs(cur-prev) / average(cur+prev)
max_err = zeros(size(out,1),1);
for qq = 1:size(out,1)
    this = x2(qq,:);
    last = x1(qq,:);
    max_err(qq) = log10( max(abs(this-last)) ./ ...
         (trapz(in2, abs(0.5*(this+last)) /( lims(2)-lims(1) ) )) );
end

%% Plotting Images, if Desired
if isempty(images)
    return
end
for nn = 1:length(images)
    plotflag = 1;
    fignum = images(nn);
    figure(fignum)
    scatter(in2,x1(nn,:),2,'ok','filled')
    hold on
    scatter(in2,x2(nn,:),2,'or','filled')
    hold off
    set(gcf,'color','w')
    xlabel('time / s')
    ylabel(['var(',num2str(nn),')'])
    legend('lastrun','currun')
end

% Make difference plots for all previously plotted variables
if plotflag ==1
    figure(fignum+1)
    hold on
    for pp = 1:size(out,1)
        leg = ['var(',num2str(pp),')'];
        scatter(inn2,xx2(pp,:)-xx1(pp,:),2,'o','DisplayName',leg)
    end
    set(gcf,'color','w')
    xlabel('time / s')
    ylabel('\Delta')
    legend show
    hold off
end

end

