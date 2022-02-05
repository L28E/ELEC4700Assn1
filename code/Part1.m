%=======================================================================
% ELEC 4700 Assignment 1
% Nicholas Ramkhalawansingh

% Part 1
%=======================================================================
clear
close all

m_0=9.10938e-31;        % electron rest mass (kg)
m_n=0.26*m_0;           % electron effective mass (kg)
T=300;                  % Temperature (K)
k_b=1.380649e-23;       % Boltzmann Constant (J/K)

% a
V_th=sqrt(2*k_b*T/m_n);   % Thermal velocity (m/s)
disp("Thermal Velocity: " + V_th + " m/s")

% b
tau_mn=0.2e-12; % Mean time between collisions 
l=V_th*tau_mn;  % Mean free path
disp("Mean Free Path: " + l + " m")

% c
num_electrons=1000;
num_steps=1000;
num_traces=5;
ymax=100e-9;
xmax=200e-9;

% Set a time step based on a 1/100th spacial step
%(min(xlim,ylim)*0.01/V_th) = 7.5619e-15 s, we'll round to 8 femtoseconds
%dt=8e-15;
dt=4e-15;

% Generate random electron positions
Px=rand(1,num_electrons).*xmax;
Py=rand(1,num_electrons).*ymax;

% Generate random electron directions (velocity is constant for now)
phi=rand(1,num_electrons)*2*pi;
Vx=V_th*cos(phi);
Vy=V_th*sin(phi);

% Randomly select some electrons to follow
tracked_indices=randperm(num_electrons,num_traces);

% Make vectors to store the paths of those electons, and temperature
X=zeros(num_traces,num_steps);
Y=zeros(num_traces,num_steps);
t=zeros(1,num_steps);
temp=zeros(1,num_steps);

figure()

for k=2:num_steps
    % Update positions
    Px=Px+Vx*dt;
    Py=Py+Vy*dt;
    
    % Electrons leaving lateral bounds come back in to preserve density
    Px(Px<0)=xmax+Px(Px<0);
    Px(Px>xmax)=Px(Px>xmax)-xmax;
    
    % Electrons reflect off upper and lower bounds
    beyond_upper=Py>ymax;
    beyond_lower=Py<0;
    Vy(beyond_lower|beyond_upper)=-Vy(beyond_lower|beyond_upper);
    Py(beyond_lower)=-Py(beyond_lower);
    Py(beyond_upper)=-Py(beyond_upper)+2*ymax;    
    
    % Update the tracked electrons and temperature 
    t(k)=t(k-1)+dt;
    X(:,k)=Px(tracked_indices);
    Y(:,k)=Py(tracked_indices);    
    temp(k)=(sum(Vx(:).^2)+sum(Vy(:).^2))*m_n/k_b/2/num_electrons;
    
    % Plot electron trajectories        
    plot(X(1,1:k),Y(1,1:k),".",X(2,1:k),Y(2,1:k),".",X(3,1:k),Y(3,1:k),".",...
        X(4,1:k),Y(4,1:k),".",X(5,1:k),Y(5,1:k),".",'MarkerSize',4) % hide the discontinuities by using dots :)
    title("Calculated Temperature: " + temp(k))   
   
    pause(0.00001)    
end

% Plot temperature
figure()
plot(t,temp)
axis([0 dt*num_steps 0 400])
title("Temperature vs. Time")
ylabel("Temperature (K)")
xlabel("Time (s)")