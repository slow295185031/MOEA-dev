function test1()
%局部搜素寻找二维四分之一圆平衡态
    N=100; %点数量
    generation = 200;
    x=rand(N,1)*pi/2;%根据随机求面均匀分布，先生成一个初始状态
    r0=[cos(x),sin(x)];
    v0=zeros(size(r0));
    G=1e-5;%斥力常数，试验这个值比较不错
    
    h=plot(r0(:,1),r0(:,2),'o');hold on;%画结果
    axis equal;
    axis([0 1 0 1]);

    for gen = 1 : generation
        [rn,vn]=countnext(r0,v0,G);%更新状态
        r0=rn;v0=vn;
        set(h,'xdata',rn(:,1),'ydata',rn(:,2));
        pause();
        drawnow;
    end
end

function [rn vn]=countnext(r,v,G) %更新状态的函数
%r存放每点的x，y，z数据，v存放每点的速度数据
num=size(r,1);
dd=zeros(2,num,num); %各点间的矢量差
for m=1:num-1
    for n=m+1:num
        dd(:,m,n)=(r(m,:)-r(n,:))';
        dd(:,n,m)=-dd(:,m,n);
    end
end
L=sqrt(sum(dd.^2,1));%各点间的距离
L(L<1e-2)=1e-2; %距离过小限定
F=sum(dd./repmat(L.^3,[2 1 1]),3)';%计算合力

Fr=r.*repmat(dot(F,r,2),[1 2]); %计算合力径向分量
Fv=F-Fr; %切向分量

rn=r+v;  %更新坐标
rn(rn<0) = 0;
rn(rn>1) = 1;
rn=rn./repmat(sqrt(sum(rn.^2,2)),[1 2]);
vn=v+G*Fv;%跟新速度
end