function [P,V]=init_UAV_ROI(R,v_min,v_max,N)
V=zeros(N,2);
if v_min<0.2
    v_min=0.2;
end
Vmag = v_min + (v_max - v_min)*rand(N,1);
Vmag(Vmag<v_min)=v_min;
Vmag=v_max;
theta=2*pi*rand(N,1);
u=rand(N,1);

x=R.*cos(theta);
y=R.*sin(theta);

% max_phi=max(phi)
% min_phi=min(phi)
P=[x,y];

%--- 速度方向：全空间均匀，再强制朝原点 ---
theta_v = 2*pi*rand(N,1);

Vdir = [cos(theta_v), sin(theta_v)];

% 若点积>0（朝外），翻转方向即可保证朝内
outward = sum(P.*Vdir,2) > 0;
Vdir(outward,:) = -Vdir(outward,:);

V=Vdir.*Vmag;
end