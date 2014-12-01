function MAIN()
% This is for testing the distribution
N=200; %点数量
a=rand(N,1)*pi/2;%根据随机求面均匀分布，先生成一个初始状态
b=asin(rand(N,1));
r0=[cos(a).*cos(b),sin(a).*cos(b),sin(b)];
v0=zeros(size(r0));
G=1e-5;%斥力常数，试验这个值比较不错

h=plot3(r0(:,1),r0(:,2),r0(:,3),'o');hold on;%画结果
axis equal;
axis([0 1 0 1 0 1]);

for ii=1:1000%模拟200步，一般已经收敛，加点判断条件其实可以在之前退出
    [rn,vn]=countnext(r0,v0,G);%更新状态
    r0=rn;v0=vn;
    set(h,'xdata',rn(:,1),'ydata',rn(:,2),'zdata',rn(:,3));
    pause();
    drawnow;
end
% dt = DelaunayTri(rn);  %利用Delaunay将点划分为空间4面体
% [ch] = convexHull(dt); %利用convexHull求凸包表面和凸包体积
% trisurf(ch,rn(:,1),rn(:,2),rn(:,3),'FaceColor','c');%画凸多面体网格
hold off;
end

function [rn vn]=countnext(r,v,G) %更新状态的函数
%r存放每点的x，y，z数据，v存放每点的速度数据
num=size(r,1);
dd=zeros(3,num,num); %各点间的矢量差
for m=1:num-1
    for n=m+1:num
        dd(:,m,n)=(r(m,:)-r(n,:))';
        dd(:,n,m)=-dd(:,m,n);
    end
end
L=sqrt(sum(dd.^2,1));%各点间的距离
L(L<1e-2)=1e-2; %距离过小限定
F=sum(dd./repmat(L.^3,[3 1 1]),3)';%计算合力

Fr=r.*repmat(dot(F,r,2),[1 3]); %计算合力径向分量
Fv=F-Fr; %切向分量

rn=r+v;  %更新坐标
rn(rn<0) = 0;
rn(rn>1) = 1;
rn=rn./repmat(sqrt(sum(rn.^2,2)),[1 3]);
vn=v+G*Fv;%跟新速度
end