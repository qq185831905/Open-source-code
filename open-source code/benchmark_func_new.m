function f=benchmark_func_new(x,func_num)
if func_num==1
    fhd=str2func('func01');  %[-100, 100]
elseif func_num==2
    fhd=str2func('func02'); %[-10;10]
elseif func_num==3 
    fhd=str2func('func03'); %[-100, 100]
elseif func_num==4
    fhd=str2func('func04');%[-100, 100]
elseif func_num==5
    fhd=str2func('func05');%[-30, 30]
elseif func_num==6
    fhd=str2func('func06'); %[-100, 100]
elseif func_num==7
    fhd=str2func('func07'); %[-1.28, 1.28]
elseif func_num==8
    fhd=str2func('func08'); %[-500, 500]
elseif func_num==9
    fhd=str2func('func09'); %[-5.12, 5.12]
elseif func_num==10
    fhd=str2func('func10'); %[-32, 32]
elseif func_num==11
    fhd=str2func('func11'); %[-600, 600]
elseif func_num==12
    fhd=str2func('func12'); %[-50, 50]
elseif func_num==13
    fhd=str2func('func13'); %[-50, 50]
elseif func_num==14
    fhd=str2func('func14'); %[-65.536;65.536]
elseif func_num==15
    fhd=str2func('func15'); %[-5, 5]
elseif func_num==16
    fhd=str2func('func16'); %[-5;5]
elseif func_num==17
    fhd=str2func('func17'); %[-5 10,0 15]
elseif func_num==18
    fhd=str2func('func18'); %[-2;2]
elseif func_num==19
    fhd=str2func('func19'); %[0, 1]
elseif func_num==20
    fhd=str2func('func20'); %[0;1]
elseif func_num==21
    fhd=str2func('func21'); %[0, 10]
elseif func_num==22
    fhd=str2func('func22'); %[0, 10]
elseif func_num==23
    fhd=str2func('func23'); %[0, 10]
end
f=feval(fhd,x);
%%%
function fit=func01(x)
% fit=sum(x.*x,2);
for i=1:size(x,1)
    s=0;
    for j=1:size(x,2)
        s=s+x(i,j)^2;
    end
    fit(i,1)=s;
end

%%%
function fit=func02(x)
x=abs(x);
fit=sum(x,2)+prod(x,2);
%%%
function fit=func03(x)
[popsize,D]=size(x);
fit=[];
for k=1:popsize
    ft=0;
    tmp=0;
    for i=1:D
        tmp=tmp+x(k,i);
        ft=ft+tmp*tmp;
    end
    fit=[fit;ft];
end
%%%
function fit=func04(x)
fit=max(abs(x),[],2);
%%%
function fit=func05(x)
[popsize,D]=size(x);
fit=[];
for k=1:popsize
    xk=x(k,:);
    yk=xk(1:D-1);
    zk=xk(2:D);
    fk=100*sum((zk - yk.^2).^2)+sum((yk-1).^2);
    fit=[fit;fk];
end
%%%
function fit=func06(x)
x=floor(x+0.5);
fit=sum(x.*x,2);
%%%
function fit=func07(x)
[popsize,D]=size(x);
x=x.^4;
index=ones(popsize,1)*(1:D);
fit=sum(index.*x,2);
fit=fit+rand(popsize,1);
%%%
function fit=func08(x)
[popsize,D]=size(x);
xt=sin(sqrt(abs(x)));
fit=418.9829*D-sum(x.*xt,2);
%%%
function fit=func09(x)
ft=x.*x-10*cos(2*pi*x)+10;
fit=sum(ft,2);
%%%
function fit=func10(x)
t1=-0.2*sqrt(mean(x.*x,2));
t2=mean(cos(2*pi*x),2);
fit=-20*exp(t1)-exp(t2)+20+exp(1);
%%%
function fit=func11(x)
[popsize,D]=size(x);
index=ones(popsize,1)*sqrt(1:D);
y=cos(x./index);
fit=sum(x.*x,2)/4000-prod(y, 2)+1;
%%%
function f = func12(x)
[popsize,D]=size(x);
f=[];
y=1+(x+1)/4;
for k=1:popsize
    xk=x(k,:);
    yk=y(k,:);
    t1=0;
    for i=1:(D-1)
        t1=t1+(yk(i)-1)^2*(1+10*(sin(pi*yk(i+1)))^2);
    end
    fk=(10*(sin(pi*yk(1)))^2+t1+(yk(D)-1)^2)*pi/D+u(xk,10,100,4);
    f=[f;fk];
end
%%%
function f=func13(x)
[popsize,D]=size(x);
f=[];
y=1+(x+1)/4;
for k=1:popsize
    xk=x(k,:);
    yk=y(k,:);
    t1=0;
    for i=1:(D-1)
        t1=t1+(xk(i)-1)^2*(1+(sin(3*pi*xk(i+1)))^2);
    end
    fk=0.1*((sin(pi*3*xk(1)))^2+t1+(xk(D)-1)^2*(1+(sin(2*pi*xk(D)))^2))+u(xk,5,100,4);
    f=[f;fk];
end
%%%
function result=u(xk,a,k,m)
xk=abs(xk);
index=find(xk>a);
result=sum(k.*(xk(index)-a).^m);
%%%
function fit=func14(x)
[popsize,D]=size(x);
a=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
fit=[];
for k=1:popsize
    xk=x(k,:);
    t1=0;
    for j=1:25
        t1=t1+1/(j+(xk(1)-a(1,j))^6+(xk(2)-a(2,j))^6);
    end
    result=1/(1/500+t1);
    fit=[fit;abs(result-0.998)];
end
%%%
function fit=func15(x)
[popsize,D]=size(x);
a=[0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b=[1/0.25 1/0.5 1/1 1/2 1/4 1/6 1/8 1/10 1/12 1/14 1/16];
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=0;
    for i=1:11
        result=result+(a(i)-(xk(1)*(b(i)^2+b(i)*xk(2)))/(b(i)^2+b(i)*xk(3)+xk(4)))^2;
    end
    fit=[fit;abs(result-0.0003075)];
end
%%%
function fit=func16(x)
[popsize,D]=size(x);
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=4*xk(1)^2-2.1*xk(1)^4+xk(1)^6/3+xk(1)*xk(2)-4*xk(2)^2+4*xk(2)^4;
    fit=[fit;abs(result+1.0316285)];
end
%%%
function fit=func17(x)
[popsize,D]=size(x);
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=(xk(2)-5.1/(4*(pi^2))*xk(1)^2+5*xk(1)/pi-6)^2+10*(1-1/(8*pi))*cos(xk(1))+10;
    fit=[fit;abs(result-0.397887)];
end
%%%
function fit=func18(x)
[popsize,D]=size(x);
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=(1+(xk(1)+xk(2)+1)^2*(19-14*xk(1)+3*xk(1)^2-14*xk(2)+6*xk(1)*xk(2)+3*xk(2)^2))*(30+(2*xk(1)-3*xk(2))^2*(18-32*xk(1)+12*xk(1)^2+48*xk(2)-36*xk(1)*xk(2)+27*xk(2)^2));
    fit=[fit;abs(result-3)];
end
%%%
function fit=func19(x)
[popsize,D]=size(x);
relfa=[1 1.2 3 3.2];
aa=[3 10 30; 0.1 10 35; 3 10 30; 0.1 10 35];
p=[0.3689 0.1170 0.2673; 0.4699 0.4387 0.747; 0.1091 0.8732 0.5547; 0.03815 0.5743 0.8828];
fit=[];
for i=1:popsize
    s=0;
    for k=1:4
        ss=0;
        for j=1:3
            ss=ss+aa(k,j)*(x(i,j)-p(k,j))^2;
        end
        s=s+relfa(k)*exp(-ss);
    end
    fit=[fit;abs(-s+3.862782)];
end
%%%
function fit=func20(x)
%  [popsize,D]=size(x);
%  relfa=[1 1.2 3 3.2];
%  bb=[10 3 17 3.5 1.7 8;0.05 10 17 0.1 8 14;3 3.5 1.7 10 17 8;17 8 0.05 10 0.1 14];
%  Q=[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;0.2348 0.1451 0.3522 0.2883 0.3047 0.665;0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
%  fit=[];
%  for i=1:popsize
%      s=0;
%      for k=1:4
%          ss=0;
%          for j=1:6
%              ss=ss+bb(k,j)*(x(i,j)-Q(k,j))^2;
%          end
%          s=s+relfa(k)*exp(-ss);
%      end
%      fit=[fit;abs(-s+3.32237)];
%  end
[popsize,D]=size(x);
relfa=[1 1.2 3 3.2];
bb=[10 3 17 3.5 1.7 8;0.05 10 17 0.1 8 14;3 3.5 1.7 10 17 8;17 8 0.05 10 0.1 14];
Q=[0.1312 0.1696 0.5569 0.0124 0.8283 0.5886;0.2329 0.4135 0.8307 0.3736 0.1004 0.9991;0.2348 0.1451 0.3522 0.2883 0.3047 0.665;0.4047 0.8828 0.8732 0.5743 0.1091 0.0381];
fit=[];
for i=1:popsize
    s=0;
    for k=1:4
        ss=0;
        for j=1:6
            ss=ss+bb(k,j)*(x(i,j)-Q(k,j))^2;
        end
        s=s+relfa(k)*exp(-ss);
    end
    fit=[fit;-s+3.32237];
end
%%%
function fit=func21(x)
[popsize,D]=size(x);
a=[4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3; 8 1 8 1; 6 2 6 2; 7 3.6 7 3.6];
c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=0;
    for i=1:5
        result=result-1/((xk-a(i,:))*((xk-a(i,:))')+c(i));
    end
    fit=[fit;abs(result+10.1532)];
end
%%%
function fit=func22(x)
[popsize,D]=size(x);
a=[4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3; 8 1 8 1; 6 2 6 2; 7 3.6 7 3.6];
c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=0;
    for i=1:7
        result=result-1/((xk-a(i,:))*((xk-a(i,:))')+c(i));
    end
    fit=[fit;abs(result+10.4029)];
end
%%%
function fit=func23(x)
[popsize,D]=size(x);
a=[4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7; 2 9 2 9; 5 5 3 3; 8 1 8 1; 6 2 6 2; 7 3.6 7 3.6];
c=[0.1 0.2 0.2 0.4 0.4 0.6 0.3 0.7 0.5 0.5];
fit=[];
for k=1:popsize
    xk=x(k,:);
    result=0;
    for i=1:10
        result=result-1/((xk-a(i,:))*((xk-a(i,:))')+c(i));
    end
    fit=[fit;abs(result+10.5364)];
end



