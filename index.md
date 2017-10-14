# 数值计算第一次上机报告

###拉格朗日插值多项式matlab实现
* #### **程序代码**
~~~matlab
function [C,L] = lagran( X,Y )
%lagran函数：
%   输入 - X是插值节点的矩阵
%        - Y是插值节点对应函数值的矩阵
%   输出 - C是拉格朗日插值多项式系数按降幂排列矩阵
%        - L是插值基函数多项式系数矩阵
w=length(X);
n=w-1;
L=zeros(w,w);
for k=1:n+1
    V=1;
    for j=1:n+1
        if k~=j
            V=conv(V,poly(X(j)))/(X(k)-X(j));
        end
    end
    L(k,:)=V;
end
C=Y*L;
~~~

*注释：*

>p = poly(r)（其中 r 是矢量）返回多项式的系数，其中多项式的根是 r 的元素。
>w = conv(u,v) 返回矢量 u 和 v 的卷积。如果 u 和 v 是多项式系数的矢量，对其卷积与将这两个多项式相乘等效。
*  #### **结果输出**

下面以课本第50面计算实习题1为例。下同。
~~~matlab
>> X=[0.2:0.2:1.0];
>> Y=[0.98 0.92 0.81 0.64 0.38];
>> C=lagran(X,Y)

C =

   -0.5208    0.8333   -1.1042    0.1917    0.9800
>> x=0.2:.01:1;y=polyval(C,x);
plot(x,y,'b',X,Y,'o')
~~~

*注释：*

>y = polyval(p,x) 返回在 x 处计算的 n 阶多项式的值。输入参数 p 是长度为 n+1 的矢量，其元素是按要计算的多项式降幂排序的系数。
>即$y = p_{1}x^n + p_{2}x^{n-1} + … + p_{n}~x + p_{n+1}$

输出图像如下：

![Lagran](http://oxtjshci8.bkt.clouddn.com/lagran.bmp)

###牛顿插值多项式matlab实现
* #### **程序代码**
~~~matlab
function [ C,D ] = newpoly( X,Y )
%newpoly函数：
%   输入 - X是插值节点的矩阵
%        - Y是插值节点对应函数值的矩阵
%   输出 - C是牛顿插值多项式系数按降幂排列矩阵
%        - D是差分表矩阵
n=length(X);
D=zeros(n,n);
D(:,1)=Y';
%开始计算差分表矩阵
for j=2:n
    for k=j:n
        D(k,j)=(D(k,j-1)-D(k-1,j-1))/(X(k)-X(k-j+1));
    end
end
%利用嵌套乘法形式求牛顿插值多项式
C=D(n,n);
for k=(n-1):-1:1
    C=conv(C,poly(X(k)));
    m=length(C);
    C(m)=C(m)+D(k,k);
end
~~~
*  #### **结果输出**

~~~matlab
>> X=[0.2:0.2:1.0];
>> Y=[0.98 0.92 0.81 0.64 0.38];
>> C=newpoly(X,Y)

C =

   -0.5208    0.8333   -1.1042    0.1917    0.9800
>> x=0.2:.01:1;y=polyval(C,x);
>> plot(x,y,'b',X,Y,'o')
~~~

输出图像如下：

![Newpoly](http://oxtjshci8.bkt.clouddn.com/Newpoly.bmp)

###三次样条插值多项式matlab实现
* #### **程序代码**
~~~matlab
function S = csfit( X,Y,dx0,dxn )
%csfit函数：
%   输入 - X是插值节点的矩阵
%        - Y是插值节点对应函数值的矩阵
%        - dx0，dxn是边界条件，若为自然白楠姐条件，则为0
%   输出 - S的每行分别是对应插值点区间内插值函数的系数
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);

B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));

for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end

M(N)=U(N-1)/B(N-1);

for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end

M(1)=3*(D(1)-dx0)/H(1)-M(2)/2;
M(N+1)=3*(dxn-D(N))/H(N)-M(N)/2;

for k=N-1:-1:0
    S(N-k,1)=(M(N-k+1)-M(N-k))/(6*H(N-k));
    S(N-k,2)=M(N-k)/2;
    S(N-k,3)=D(N-k)-H(N-k)*(2*M(N-k)+M(N-k+1))/6;
    S(N-k,4)=Y(N-k);
end
~~~

*注释：*

>这里的三次样条函数最终形式和书上不同，具体为
>$S(x)=S_{i+1,1}(x-x_{i})^3+S_{i+1,2}(x-x_{i})^2+S_{i+1,3}(x-x_{i})+S_{i+1,4},$
>$x\in[x_{i},x_{i+1}], \quad i=0,1,2...n-1$
*  #### **结果输出**

~~~matlab
>> X=[0.2:0.2:1.0];
>> Y=[0.98 0.92 0.81 0.64 0.38];
>> dx0=0;dxn=0;
>> S = csfit( X,Y,dx0,dxn )

S =

    2.5446   -2.0089         0    0.9800
    1.1161   -0.4821   -0.4982    0.9200
   -8.2589    0.1875   -0.5571    0.8100
   28.1696   -4.7679   -1.4732    0.6400
>> for i=1:4
	xi=(0.2+0.2*(i-1)):0.01:(0.2+0.2*i);
	yi=polyval(S(i,:),xi-X(i));
	plot(xi,yi,'b')
	hold on
end
>> plot(X,Y,'o')
~~~

输出图像如下：

![Csfit](http://oxtjshci8.bkt.clouddn.com/Csfit.bmp)

###计算实习题1-图像汇总
  由于拉格朗日插值多项式与牛顿插值多项式本质上相同，故略去牛顿多项式的作图
~~~matlab
>> X=[0.2:0.2:1.0];
>> Y=[0.98 0.92 0.81 0.64 0.38];
>> dx0=0;dxn=0;
>> S = csfit( X,Y,dx0,dxn );
>> C=lagran(X,Y);
>> for i=1:4
	xi=(0.2+0.2*(i-1)):0.01:(0.2+0.2*i);
	yi=polyval(S(i,:),xi-X(i));
	plot(xi,yi,'b')
	hold on
end
>> x5=0.2:.01:1;y5=polyval(C,x5);
>> plot(x5,y5,'r',X,Y,'o')
~~~

输出图像如下，其中蓝线代表三次样条插值函数，红线代表拉格朗日插值函数

![Lagran+Csfit](http://oxtjshci8.bkt.clouddn.com/Lagran+Csfit.bmp)


###计算实习题2-图像汇总
*  #### **n=10**
~~~matlab
syms x;
runge(x)=1/(1+25*x^2);
X=[-1:0.2:1];
Y=zeros(1,11);
for k=1:11
	Y(k)=runge(X(k));
end
dx0=diff(runge(-1));
dxn=diff(runge(1));

L=lagran(X,Y);
S=csfit(X,Y,dx0,dxn);

for i=1:10
	xi=(-1+0.2*(i-1)):0.01:(-1+0.2*i);
	yi=polyval(S(i,:),xi-X(i));
	plot(xi,yi,'b')
	hold on
end

x11=-1:0.01:1;y11=polyval(L,x11);
plot(x11,y11,'g')
hold on

x12=-1:0.01:1;y12=runge(x12);
plot(x12,y12,'-r')
hold on

plot(X,Y,'or')
~~~

图像如下，其中绿线代表拉格朗日插值函数，红线代表龙格函数，蓝线代表三次样条插值函数
![n=10](http://oxtjshci8.bkt.clouddn.com/n=10.bmp)

*  #### **n=20**
~~~matlab
syms x;
runge(x)=1/(1+25*x^2);
X=[-1:0.1:1];
Y=zeros(1,21);
for k=1:21
	Y(k)=runge(X(k));
end
dx0=diff(runge(-1));
dxn=diff(runge(1));

L=lagran(X,Y);
S=csfit(X,Y,dx0,dxn);

for i=1:20
	xi=(-1+0.1*(i-1)):0.01:(-1+0.1*i);
	yi=polyval(S(i,:),xi-X(i));
	plot(xi,yi,'b')
	hold on
end

x21=-1:0.01:1;y21=polyval(L,x21);
plot(x21,y21,'g')
hold on

x22=-1:0.01:1;y22=runge(x22);
plot(x22,y22,'-r')
hold on

plot(X,Y,'or')
~~~
图像如下，其中绿线代表拉格朗日插值函数，红线代表龙格函数，蓝线代表三次样条插值函数
![n=20](http://oxtjshci8.bkt.clouddn.com/n=20.bmp)

可以看出，随着次数的升高，拉格朗日插值函数的病态性质越明显

---

**林嘉恒 15级信息与计算科学 2015111429**

>*主要参考书目：*
>John M, Kurtis F 数值方法（MATLAB版）（第四版）北京：电子工业出版社，2010

