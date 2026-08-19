clc;clear all;close all;
h_max_list=10:3000;
% fai_cover=0:pi/100:pi/2;
fai_cover=linspace(pi/2,0,180);
omiga_list=[2*pi,2*pi];%速度，rad/s
fasan_list=[10e-3 45e-3 150e-3];%发散全角
fasan_list=1e-3:1e-3:200e-3;
x_tklb_list={"10","50","90","130"};
% y_tklb_list={'0','\pi/6','\pi/3','\pi/2'};
y_tklb_list={'\pi/2','\pi/3','\pi/6','0'};
z=3000;
R=z;
t1=zeros(length(fai_cover),length(fasan_list));
omiga=omiga_list(1);
ColorList=["#808080" "#FF6B6B" "#FFC119" "#039BE5" "#5E35B1"];%深色

cmp="turbo";
linear_or_log="log";
FONT_SIZE=12;


for i=1:length(fai_cover)
    for j=1:length(fasan_list)
    D=2*z.*tan(fasan_list(j)/2);
    delta_s=D;
    k=delta_s./2./pi./R;
    func=@(x) sqrt(k.^2+sin(k.*x).^2);
    % faimin=pi/2-asin(h_max_list(j)./R)+asin(D/2./R);
    faimin=pi/2-fai_cover(i)+asin(D/2./R);
    faimax=pi/2-asin(D/2./R);
    minlim=faimin./k;
    maxlim=faimax./k;
    % minlim=(pi/2-asin(h_max_list(j)./R)+asin(D/2./R))./k;
    % maxlim=(pi/2-asin(D/2./R))./k;
    fai1=minlim*k;
    fai2=maxlim*k;
    integ=integral(func,minlim,maxlim);
    t1(i,j)=R.*integ./z./omiga+2*pi*sin(faimax)./omiga+2*pi*sin(faimin)./omiga;
    end
end
ax1=gca;
lwid=2;
% h=heatmap(t);
imagesc(fasan_list,fai_cover,t1);
colormap(cmp);
xlabel('\it{w}\rm_D (mrad)','FontSize',FONT_SIZE);
xticks([50 100 150 200]*1e-3);
xticklabels([50 100 150 200]);
yticks(linspace(0,pi/2,4));
yticklabels(y_tklb_list);
ylabel('\phi (rad)','FontSize',FONT_SIZE);
set(gca,'fontsize',FONT_SIZE);
set(gca, 'ColorScale',linear_or_log);



%%


clm = [0.5 200];     %两图共用的clim范围
clim(ax1, clm);
cb = colorbar(ax1);          % 以第二个轴为“宿主”创建colorbar
cb.Label.String = 'Scan Time (s)';
set(cb,'Ticks',[1,10,100],'FontSize',FONT_SIZE);


hold(ax1,'on');
T_iso = 2;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改
contour(ax1, fasan_list, fai_cover, t1, [T_iso T_iso],'w--', 'LineWidth', 2);
T_iso = 5;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改
contour(ax1, fasan_list, fai_cover, t1, [T_iso T_iso],'w-', 'LineWidth', 2);
T_iso = 10;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改
contour(ax1, fasan_list, fai_cover, t1, [T_iso T_iso],'w:', 'LineWidth', 2);

hold(ax1,'off');


set(gca,'fontsize',FONT_SIZE);
set(gcf,'position',[100 100 580 330])
exportgraphics(gcf, 'D:\lzx\vortex\多种方式对比小论文吗\原始图片\EA扫描时间连续发散变量.png', 'Resolution', 300);  % 设置DPI为300

