%=======================================================================
% ELEC 4700 Assignment 1
% Nicholas Ramkhalawansingh

% Part 2
%=======================================================================
clear
close all

m_0=9.10938e-31;        % electron rest mass (kg)
m_n=0.26*m_0;           % electron effective mass (kg)
T=300;                  % Temperature (K)
k_b=1.380649e-23;       % Boltzmann Constant (J/K)

V_th=sqrt(2*k_b*T/m_n);   % Thermal velocity (m/s)
tau_mn=0.2e-12; % Mean time between collisions 

num_electrons=1000;
num_steps=1000;
num_traces=500; % Increase the number of tracked electrons for MFP accuracy
ymax=100e-9;
xmax=200e-9;
dt=4e-15;
P_scat=1-exp(-dt/tau_mn);

% Generate random electron positions
Px=rand(1,num_electrons).*xmax;
Py=rand(1,num_electrons).*ymax;

% Generate random electron velocities (Normal distribution for each component of velocity)
%twiddle=1.25;  
twiddle=1; 
Vx=randn(1,num_electrons)*sqrt(k_b*T/m_n)*twiddle;
Vy=randn(1,num_electrons)*sqrt(k_b*T/m_n)*twiddle;

% Randomly select some electrons to follow
tracked_indices=randperm(num_electrons,num_traces);

% Make vectors to store the paths of those electons, and temperature
X=zeros(num_traces,num_steps);
Y=zeros(num_traces,num_steps);
t=zeros(1,num_steps);
temp=zeros(1,num_steps);

% 2D array to track the timesteps where each electron has a collision
collisions=zeros(num_electrons,num_steps);

figure()

for k=2:num_steps
    % Update positions
    Px=Px+Vx*dt;
    Py=Py+Vy*dt;
    
    % Scatter electrons
    scat=rand(1,num_electrons)<P_scat;
    Vx(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n)*twiddle;
    Vy(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n)*twiddle;    
    
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
    
    % Record the time steps where the electrons scatter
    %collisions(scat|beyond_upper|beyond_lower,k)=1;
    collisions(scat,k)=1;
    
    % Plot electron trajectories        
    plot(X(1,1:k),Y(1,1:k),".",X(2,1:k),Y(2,1:k),".",X(3,1:k),Y(3,1:k),".",...
        X(4,1:k),Y(4,1:k),".",X(5,1:k),Y(5,1:k),".",'MarkerSize',4) % hide the discontinuities by using dots :)
    title("Calculated Temperature: " + temp(k))   
   
    pause(0.00001)    
end

% Plot temperature
figure()
plot(t,temp)
axis([0 dt*num_steps 0 500])
title("Temperature vs. Time")
ylabel("Temperature (K)")
xlabel("Time (s)")

% Plot electron velocity histogram
V=sqrt(Vx(:).^2+Vx(:).^2);
figure()
histogram(V)
title('# Electrons at a Given Velocity')
xlabel('Electron Velocity (m/s)')
ylabel('# Electrons')

disp("Thermal Velocity: "+V_th)
disp("Mean Electron Velocity: "+mean(V))

% Mean free path, average for each electron and then average them all 
tot=0;
for j=1:num_traces    
    horz_dists=diff(X(j,find(collisions(j,:))));
    vert_dists=diff(Y(j,find(collisions(j,:))));    
    tot=tot+mean(sqrt(horz_dists.^2 + vert_dists.^2));
end

disp("mean free path: "+tot/num_traces);

% Mean time between collisions, average for each electron and then average them all
tot=0;
for j=1:num_electrons
    tot=tot+mean(diff(find(collisions(j,:))))*dt;
end

disp("mean time between collisions: "+tot/num_electrons);
