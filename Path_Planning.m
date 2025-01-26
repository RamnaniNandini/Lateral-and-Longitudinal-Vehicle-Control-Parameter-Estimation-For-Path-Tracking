%%
clc
clear all
%% Speed input and veh. param
%vx=5;
ds=2.5;      % Look ahead distance
%time_step=ds/vx;    % time-step
importdata input.xlsx;
speed = ans.speed;
time = ans.time;
m = 1295;
l=2.46;
lf = 0.98;
lr = l-lf;
Iz = m*lf*lr;
cf = 155494.663;
cr = 155494.663;
t1(1,1) = 0;
%% Act_dist and 1st input
lat_act=ans.lat(:,2);
long_act=ans.long(:,2);n1=size(lat_act);
z=zeros(1,1);origin = [lat_act(1), long_act(1), z];[x1,y1] = latlon2local(lat_act,long_act,z,origin);

for i=1:n1-2
K_curv(i+1,1)=2*abs((x1(i+1)-x1(i)).*(y1(i+2)-y1(i))-(x1(i+2)-x1(i)).*(y1(i+1)-y1(i))) ./sqrt(((x1(i+1)-x1(i)).^2+(y1(i+1)-y1(i)).^2).*((x1(i+2)-x1(i)).^2+(y1(i+2)-y1(i)).^2).*((x1(i+2)-x1(i+1)).^2+(y1(i+2)-y1(i+1)).^2));
end
K_curv(end+1,1)=K_curv(end,1);
% K_curv(1:253,1)=1/520;

% Haversine formula
for i_a=1:n1-1
R_Earth_act=6371;
diff_lat_act=deg2rad(lat_act(i_a+1)-lat_act(i_a));
diff_long_act=deg2rad(long_act(i_a+1)-long_act(i_a));
hav1_act=(((sin(diff_lat_act/2)))^2)+(cos(deg2rad(lat_act(i_a+1)))*cos(deg2rad(lat_act(i_a)))*((sin(diff_long_act/2))^2));
hav2_act=2*R_Earth_act*atan2(sqrt(hav1_act),sqrt(1-hav1_act));
hav_act(1)=0;
hav_act(i_a+1)=hav2_act;
end
hav_act=hav_act';
hav_act=hav_act*1000;
hav_act_cum=cumsum(hav_act);
act_dist=sum(hav_act);


curv=interp1(hav_act_cum,K_curv,ds);

%r_s1=curv*vx;    
del1(1,1) = 0;

%%
I_hdg_vbox=19.9;
I_hdg=(5*pi/2)-(I_hdg_vbox*pi/180);
hdg_vbox=ans.hdg(:,:);
I_hdg_vbox2=interp1(hav_act_cum,hdg_vbox(:,2),ds);
I_hdg_2=(5*pi/2)-(I_hdg_vbox2*pi/180);
x0 = [0;0;0;0];
curv=interp1(hav_act_cum,K_curv,ds);

K_curv_temp = curv;
k = 1;
sim_dist=0;
k=1;
sim_dist = 2.5;
distance(1,1)=sim_dist;
K1 = [0,0,0,0];
%% 
while (sim_dist<=act_dist)
vx = interp1(hav_act_cum,speed,sim_dist);
%vx = speed(k,1);
A = [0, 1, 0, 0,;0, -(cf+cr)/(m*vx), (cf+cr)/m, ((lr*cr-lf*cf)/(m*vx)); 0, 0, 0, 1; 0, (lr*cr-lf*cf)/(Iz*vx), (lf*cf-lr*cr)/Iz, -(lf^2*cf+lr^2*cr)/(Iz*vx)];
B = [0; cf/m; 0; lf*cf/Iz];
C=eye(size(A));
% C = [1 0 0 0;0 0 0 0; 0 0 1 0;0 0 0 0];
D = zeros(size(B));
sys = ss(A, B, C, D);

% Define Q and R matrices
Q = [200,0,0,0;
    0,0,0,0;
    0,0,10000,0;
    0,0,0,0];
R = 1;
% Get optimal gain
[K,S,P] = lqr(sys,Q,R);
K1 = [K1;K];
%Update state-space matrices based on optimal gain
A1 = A - B*K;
B2 = [0; ((lr*cr-lf*cf)/(m*vx))-vx; 0;-(lf^2*cf+lr^2*cr)/(Iz*vx) ];
time_step = ds/vx;
t1(k+1,1) = time_step;
tspan = [0 time_step];
A=A1;
u=K_curv_temp*vx;
%x0 = [0; 0; 0; 0];
[t,y] = ode45(@(t,y) CDE2(t,y,A,B2,u), tspan, x0);

% updating variables
ecg_temp=y(end,1);
e_dot_cg_temp=y(end,2);
theta_e_temp=y(end,3);
theta_dot_e_temp=y(end,4);
%vy_temp=e_dot_cg_temp-vx*theta_e_temp;
%r_temp=theta_dot_e_temp+u;   % check this once


% final variables
ecg(k+1,1)=ecg_temp;
e_dot_cg(k+1,1)=e_dot_cg_temp;
theta_e(k+1,1)=theta_e_temp;
theta_dot_e(k+1,1)=theta_dot_e_temp;
%vy(k+1,1)=vy_temp;
%r(k+1,1)=r_temp;

%time_temp=[0; time_step];
%hdg_temp=cumtrapz(time_temp,[r(k,1);r(k+1,1)]);
%hdg_temp=(-1)*hdg_temp;
%psi_temp=psi(k,1)+hdg_temp(end,1)+atan(vy_temp/vx);   % heading with respect to east in counter clockwise direction
%psi(k+1,1)=psi_temp;

%R = [(cos(theta)*cos(psi(k+1))),(sin(phi)*sin(theta)*cos(psi(k+1)))-(cos(phi)*sin(psi(k+1)));(cos(theta)*sin(psi(k+1))),(sin(phi)*sin(theta)*sin(psi(k+1)))+(cos(phi)*cos(psi(k+1)))];
%Vx_Vy=[vx;vy_temp];
%G_Vx_Vy=R*Vx_Vy;
%G_Vx(k+1,1)=G_Vx_Vy(1,1);
%G_Vy(k+1,1)=G_Vx_Vy(2,1);
%X_temp=cumtrapz(time_temp,[G_Vx(k,1);G_Vx(k+1,1)]);
%Y_temp=cumtrapz(time_temp,[G_Vy(k,1);G_Vy(k+1,1)]);
%X(k+1,1)=X(k)+X_temp(end,1);
%Y(k+1,1)=Y(k)+Y_temp(end,1);
%sim_dist_temp=((X(k+1)-X(k))^2+(Y(k+1)-Y(k))^2)^0.5;
%sim_dist_1(k+1,1)=sim_dist_temp;
%sim_dist=sum(sim_dist_1);

K_curv_temp=interp1(hav_act_cum,K_curv,sim_dist+ds);
sim_dist = sim_dist+ds;
distance = [distance;sim_dist];

x0=[ecg_temp;e_dot_cg_temp;theta_e_temp;theta_dot_e_temp];
del = -K*x0;
del1 = [del1;del];
k=k+1;disp(k)
end
%% 
t = t1;
t_end = time_step*(k-1);
time=[0:time_step:t_end];time=time';
%distance=cumsum(sim_dist_1);
X_int=interp1(hav_act_cum,x1,distance);
Y_int=interp1(hav_act_cum,y1,distance);
K_curv_int=interp1(hav_act_cum,K_curv,distance);
r_des=vx*K_curv_int;
psi_des=cumtrapz(time,r_des);

%%
I_hdg_vbox=19.9;
I_hdg=(5*pi/2)-(I_hdg_vbox*pi/180);I_hdg(1:258)=I_hdg;I_hdg=I_hdg';
psi_des_new=I_hdg-psi_des;
%%
X_sim=X_int-ecg.*sin(theta_e+psi_des_new);
Y_sim=Y_int+ecg.*cos(theta_e+psi_des_new);

%%
plot(x1,y1);hold on;
plot(X_sim,Y_sim);

%%
x=[ecg,e_dot_cg,theta_e,theta_dot_e];
delta=(-K*x');
figure;
plot(time,r_des);
hold on;plot(time,delta);

%% 

z=zeros(1,1);
origin = [lat_act(1), long_act(1), z];
[lat_mod,long_mod] = local2latlon(X_sim,Y_sim,z,origin);

hold off
figure;
geoplot(lat_mod,long_mod,'g.')
hold on

% Plotting actual lat-long
geoplot(lat_act,long_act,'r.')
geobasemap satellite
%% Comparison bw del_sim and del_act
del_act=(ans.Sheet1(:,2)*pi/180)/12;
del_act2=(ans.Sheet1(:,3)*pi/180)/12;
time_act=ans.Sheet1(:,1);
figure;

plot(time,del1);hold on;
plot(time_act,del_act2);
%% 
function dydt=CDE2(t1,y,A,B2,u)
dydt=zeros(4,1);

dydt(1) = y(2);
dydt(3)=y(4);
y(2) = A(1,1)*y(1)+A(1,2)*y(2)+A(1,3)*y(3)+A(1,4)*y(4)+B2(1,1)*u;
dydt(2) = A(2,1)*y(1)+A(2,2)*y(2)+A(2,3)*y(3)+A(2,4)*y(4)+B2(2,1)*u;
y(4) = A(3,1)*y(1)+A(3,2)*y(2)+A(3,3)*y(3)+A(3,4)*y(4)+B2(3,1)*u;
dydt(4) = A(4,1)*y(1)+A(4,2)*y(2)+A(4,3)*y(3)+A(4,4)*y(4)+B2(4,1)*u;
end
