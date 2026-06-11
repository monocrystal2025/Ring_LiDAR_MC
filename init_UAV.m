% [PP,VV]=init_UAV1(100,0.5,5,10000);
function [P,V]=init_UAV(R,v_min,v_max,N)
V=zeros(N,3);
if v_min<0.2
    v_min=0.2;
end
Vmag = v_min + (v_max - v_min)*rand(N,1);
Vmag(Vmag<v_min)=v_min;
Vmag=v_max;
theta=2*pi*rand(N,1);
u=rand(N,1);
phi=acos(u);
phi1=phi;
x=R*sin(phi).*cos(theta);
y=R*sin(phi).*sin(theta);
z=R*cos(phi);
% max_phi=max(phi)
% min_phi=min(phi)
P=[x,y,z];

%--- 速度方向：全空间均匀，再强制朝原点 ---
theta_v = 2*pi*rand(N,1);
u_v = 2*rand(N,1)-1;        % cos(phi_v) ~ U(-1,1)
phi_v = acos(u_v);          % 弧度
Vdir = [sin(phi_v).*cos(theta_v), sin(phi_v).*sin(theta_v), cos(phi_v)];

% 若点积>0（朝外），翻转方向即可保证朝内
outward = sum(P.*Vdir,2) > 0;
Vdir(outward,:) = -Vdir(outward,:);

V=Vdir.*Vmag;
end