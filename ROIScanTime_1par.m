clear all;close all;clc;
z=50:3000;
z=1000;
% a=[100 500];
w=5:0.2:45;
a=1000*tand(w);
tao=1/5e3;     %50kHz重频激光下的脉冲间隔时间
omiga=[2*pi 2*pi];
ColorList=["#979998" "#315A89" "#0787C3" "#B42B22" "#EC3232"];
ColorList=["#808080" "#00897B" "#0288D1" "#512DA8" "#D81B60" "#E07000" "#B88300"];
% fasan_list=[10e-3 15e-3 45e-3 127e-3 15e-3 45e-3 127e-3];%发散全角
fasan_list=1e-3:0.2e-3:200e-3;
D=2*z.*tan(fasan_list./2);%不同距离下光场直径
L1=zeros(length(w),length(fasan_list));
L2=L1;


%fasan_guangban_list=30e-3;%光斑发散角

for i=1:length(a)
    for j=1:length(fasan_list)
        shangxian=2*pi*(a(i)-0.5.*D(j))./D(j);
        shangxian=2*pi*(a(i)+0.5.*D(j))./D(j);
        jifen=shangxian.*sqrt(shangxian.^2+1)+log(abs(shangxian+sqrt(shangxian.^2+1)));
        L1(i,j)=D(j).*jifen./4./pi+2*pi*(a(i)-0.5.*D(j));
        L1(i,j)=D(j).*jifen./4./pi;
        if a(i)<=D(j)
            L1(i,j)=pi*D(j);
        end
        L2(i,j)=2*a(i).*ceil(2*a(i)./D(j))+D(j).*(ceil(2*a(i)./D(j))-1);
        if 2*a(i)<=D(j)
            L2(i,j)=2*a(i);
        end
    end
end
t1=L1./omiga(1)./z;
t2=L1./omiga(2)./z;
t3=L2./omiga(1)./z;
t4=L2./omiga(2)./z;
max1=max(max(t1));min1=min(min(t1));
max2=max(max(t2));min2=min(min(t2));
max3=max(max(t3));min3=min(min(t3));
max4=max(max(t4));min4=min(min(t4));
ma=max([max1,max2,max3,max4]);mi=min([min1,min2,min3,min4]);
ma=120;mi=0.08;

tl=tiledlayout(2,1);
tl.Padding="compact";
tl.TileSpacing = 'compact';
cmp="turbo";
FONT_SIZE=12;
linear_or_log='log';






ax2=nexttile(2);
imagesc(fasan_list,w*2,t2);
set(gca, 'ColorScale',linear_or_log,'FontSize',FONT_SIZE);
xlabel('\it{w}\rm_D (mrad)');
xticks([50 100 150 200]*1e-3);
xticklabels([50 100 150 200]);
yticks([10 30 50 70 90]);
ylabel('\theta (\circ)');
colormap(cmp);
clim([mi,ma]);






ax4=nexttile(1);
imagesc(fasan_list,w*2,t4);
set(gca, 'ColorScale',linear_or_log,'FontSize',FONT_SIZE);
xlabel('\it{w}\rm_D (mrad)');
xticks([50 100 150 200]*1e-3);
xticklabels([50 100 150 200]);
yticks([10 30 50 70 90]);
ylabel('\theta (\circ)');
colormap(cmp);
clim([mi,ma]);
cb=colorbar();
cb.Layout.Tile = 'east';     % 放到整个布局的右侧
cb.Label.String = 'Scan Time (s)';
set(cb,'Ticks',[0.1,1,10,100],'FontSize',FONT_SIZE);


hold(ax2,'on');hold(ax4,'on');
T_iso = 2;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改

contour(ax2, fasan_list, 2*w, t2, [T_iso T_iso],'w--', 'LineWidth', 2);

contour(ax4, fasan_list, 2*w, t4, [T_iso T_iso],'w--', 'LineWidth', 2);
T_iso = 5;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改

contour(ax2, fasan_list, 2*w, t2, [T_iso T_iso],'w-', 'LineWidth', 2);

contour(ax4, fasan_list, 2*w, t4, [T_iso T_iso],'w-', 'LineWidth', 2);
T_iso = 10;   % 例如画 t = 5 s 的等时间曲线，你可以按需要修改

contour(ax2, fasan_list, 2*w, t2, [T_iso T_iso],'w:', 'LineWidth', 2);

contour(ax4, fasan_list, 2*w, t4, [T_iso T_iso],'w:', 'LineWidth', 2);

set(gcf,"Position",[50 50 550 550]);




exportgraphics(gcf, 'D:\lzx\vortex\多种方式对比小论文吗\原始图片\ROI扫描时间_Continue.png', 'Resolution', 300);  % 设置DPI为300