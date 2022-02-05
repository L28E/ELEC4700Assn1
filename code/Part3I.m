%=======================================================================
% ELEC 4700 Assignment 1
% Nicholas Ramkhalawansingh

% Part 3
%=======================================================================
clear
close all

type=0;                 % 0 for specular, anything else for diffusive 

m_0=9.10938e-31;        % electron rest mass (kg)
m_n=0.26*m_0;           % electron effective mass (kg)
T=300;                  % Temperature (K)
k_b=1.380649e-23;       % Boltzmann Constant (J/K)
tau_mn=0.2e-12;         % Mean time between collisions 

num_electrons=1000;
num_steps=1000;
num_traces=8;
ymax=100e-9;
xmax=200e-9;
dt=4e-15;
P_scat=1-exp(-dt/tau_mn);

% Generate the boxes
box_x=0.8e-7;
box_width=0.4e-7;
box_right=box_x+box_width;
box_y=0;
box_height=0.4e-7;
box_top=box_y+box_height;
box_bottom=0.6e-7;

% Generate random electron positions
Px=rand(1,num_electrons).*xmax;
Py=rand(1,num_electrons).*ymax;

% Remove electrons from the boxes
oob=get_oob(Px,Py,box_x,box_right,box_top,box_bottom);
[M,I]=get_nearest_bound(num_electrons,Px,Py,oob,box_x,box_right,box_top,box_bottom);
[Px,Py]=remove_oob(Px,Py,oob,box_x,box_right,box_top,box_bottom,M,I);

% Generate random electron velocities (Normal distribution for each component of velocity)
Vx=randn(1,num_electrons)*sqrt(k_b*T/m_n);
Vy=randn(1,num_electrons)*sqrt(k_b*T/m_n);

% Randomly select some electrons to follow
tracked_indices=randperm(num_electrons,num_traces);

% Make vectors to store the paths of those electons, and temperature
X=zeros(num_traces,num_steps);
Y=zeros(num_traces,num_steps);
t=zeros(1,num_steps);

% 2D array to track the timesteps where each electron has a collision
collisions=zeros(num_electrons,num_steps);

figure();

for k=2:num_steps
    % Update positions
    Px=Px+Vx*dt;
    Py=Py+Vy*dt;
    
    % Scatter electrons
    scat=rand(1,num_electrons)<P_scat;
    Vx(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n);
    Vy(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n);    
    
    % Electrons leaving lateral bounds come back in to preserve density
    Px(Px<0)=xmax+Px(Px<0);
    Px(Px>xmax)=Px(Px>xmax)-xmax;
    
    % Electrons reflect off upper and lower bounds
    beyond_upper=Py>ymax;
    beyond_lower=Py<0;
    Vy(beyond_lower|beyond_upper)=-Vy(beyond_lower|beyond_upper);
    Py(beyond_lower)=-Py(beyond_lower);
    Py(beyond_upper)=-Py(beyond_upper)+2*ymax;  
    
    % Logically index electrons which are out of bounds
    oob=get_oob(Px,Py,box_x,box_right,box_top,box_bottom);    

    % Determine which bound each oob electron is closest to
    [M,I]=get_nearest_bound(num_electrons,Px,Py,oob,box_x,box_right,box_top,box_bottom);  
    
    hzn_flip=oob & M'~=0 & (I'==1|I'==2); % For oob electrons nearest the left or right edges, flip their horizontal velocity
    vert_flip=oob & M'~=0 & (I'==3|I'==4); % For oob electrons nearest the top boundary, flip their vertical velocity
    
    % Reflect & remove oob electrons from boxes
    if type==0          
        Vx(hzn_flip)=-Vx(hzn_flip);          
        Vy(vert_flip)=-Vy(vert_flip);
    else
        Vx(hzn_flip)=-sign(Vx(hzn_flip)).*abs(randn(1,nnz(hzn_flip))*sqrt(k_b*T/m_n));        
        Vy(vert_flip)=-sign(Vy(vert_flip)).*abs(randn(1,nnz(vert_flip))*sqrt(k_b*T/m_n));
    end
    [Px,Py]=remove_oob(Px,Py,oob,box_x,box_right,box_top,box_bottom,M,I);   
    
    % Update the tracked electrons and temperature 
    t(k)=t(k-1)+dt;
    X(:,k)=Px(tracked_indices);
    Y(:,k)=Py(tracked_indices);    
        
    % Record the time steps where the electrons scatter
    %collisions(scat|beyond_upper|beyond_lower,k)=1;
    collisions(scat,k)=1;
    
    % Plot electron trajectories   
    plot(X(1,1:k),Y(1,1:k),".",X(2,1:k),Y(2,1:k),".",X(3,1:k),Y(3,1:k),".",...
        X(4,1:k),Y(4,1:k),".",X(5,1:k),Y(5,1:k),".",X(6,1:k),Y(6,1:k),".",...
        X(7,1:k),Y(7,1:k),".",X(8,1:k),Y(8,1:k),".",'MarkerSize',4) % hide the discontinuities by using dots :)
    hold on;
    
    rectangle('Position',[box_x box_y box_width box_height]);
    rectangle('Position',[box_x box_bottom box_width box_height]);
    hold off;
   
    pause(0.00001)    
end

% Plot density map
figure();
binscatter(Px,Py,30)
title("Electron Density Map")
axis([0 xmax 0 ymax])
colormap('parula')

% Plot temperature
%https://www.mathworks.com/help/matlab/ref/griddata.html
figure();
temp=(Vx(:).^2+Vy(:).^2)*m_n/k_b/2;
[xq,yq] = meshgrid(0:0.05*xmax:xmax,0:0.05*ymax:ymax);
vq = griddata(Px,Py,temp',xq,yq);
s=mesh(xq,yq,vq);
s.FaceColor='interp';
view(2);
hold on;
plot3(Px,Py,temp','.r');
colorbar
title("Temperature Map")
colormap('parula')

function oob = get_oob(Px,Py,left,right,top,bottom)
    oob=(Px>=left & Px<=right & Py<=top)|(Px>=left & Px<=right & Py>=bottom);
end

function [minimum,index] = get_nearest_bound(num_electrons, Px,Py,oob,left,right,top,bottom)
    diffs=zeros(num_electrons,4);
    diffs(oob,:)=[abs(Px(oob)'-left) abs(Px(oob)'-right) abs(Py(oob)'-top) abs(Py(oob)'-bottom)];
    [minimum,index]=min(diffs,[],2);
end

function [Px,Py] = remove_oob(Px,Py,oob,left,right,top,bottom,M,I)
    Px(oob & M'~=0 & I'==1)=2*left-Px(oob & M'~=0 & I'==1);
    Px(oob & M'~=0 & I'==2)=2*right-Px(oob & M'~=0 & I'==2);
    Py(oob & M'~=0 & I'==3)=2*top-Py(oob & M'~=0 & I'==3);
    Py(oob & M'~=0 & I'==4)=2*bottom-Py(oob & M'~=0 & I'==4);
end
