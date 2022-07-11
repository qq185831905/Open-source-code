%mex cec14_func.cpp -DWINDOWS
clc
clear
func_num=1;

%总参数
pop_size=10;
iter_max=500;


% Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Function_name='F23'; 
% Load details of the selected benchmark function
[Xmin,Xmax,D,fobj]=get_functions_details30(Function_name);
Max_FES=10000*D;
        
%设定一个初始化种群
setpop=initialization(pop_size,D,Xmax,Xmin);


setpop=setpop;
MaxFes=Max_FES;
VarNumber=D;
nPop=pop_size;
T=iter_max;
LB=Xmin;
UB=Xmax;
fobj=fobj;
Costave=0;
temp=0;
%setpop 初始种群
%MaxFes 默认10000*变量维度=计算适应度的次数
%VarNumber 变量维度
%nPop 种群数量
%T 迭代最大值
%LB 下限
%UB 上限
%fobj 目标函数
%% Problem Information
VarMin = LB *ones(1,VarNumber);        % Lower bound of variable;
VarMax = UB *ones(1,VarNumber);         % Upper bound of variable;
BestCost=10000;
%% General Parameters of the Algorithm

LayerNumber = 5 ;     % Maximum number of Layers around nucleus
FotonRate = 0.1 ;     % Foton Rate for position determination of electrons

%% Counters
Iter=0;
FEs=0;

%% Initialization
Pop=[]; Cost=[];
for i=1:nPop
    % Initialize Positions
    Pop(i,:)=setpop(i,:);
    % Cost Function Evaluations
    Cost(i,1)=fobj(Pop(i,:));
    FEs=FEs+1;
end

% Sort Population
[Cost, SortOrder]=sort(Cost);
Pop=Pop(SortOrder,:); %按顺序排好
BestPop=Pop(1,:);  
MeanPop=mean(Pop); %每列的均值

% FEs=FEs1=FEs1=FEs1=FEs1;
% Iter1
% LayerNumber1
% nPop1
% Pop1
% Cost1
% FotonRate1
% temp1
% Costave1
% UB1
% LB1
% VarNumber1
% VarMin1
% VarMax1
% fobj1
% BestPop1
% MeanPop1


while FEs<MaxFes
[BestCost1,CostB1,FEs1,Iter1,LayerNumber1,nPop1,Pop1,Cost1,FotonRate1,temp1,Costave1,UB1,LB1,VarNumber1,VarMin1,VarMax1,fobj1,BestPop1,MeanPop1]=CACS(FEs1,Iter1,LayerNumber1,nPop1,Pop1,Cost1,FotonRate1,temp1,Costave1,UB1,LB1,VarNumber1,VarMin1,VarMax1,fobj1,BestPop1,MeanPop1);
[BestCost2,CostB2,FEs2,Iter2,LayerNumber2,nPop2,Pop2,Cost2,FotonRate2,temp2,Costave2,UB2,LB2,VarNumber2,VarMin2,VarMax2,fobj2,BestPop2,MeanPop2]=CACS(FEs2,Iter2,LayerNumber2,nPop2,Pop2,Cost2,FotonRate2,temp2,Costave2,UB2,LB2,VarNumber2,VarMin2,VarMax2,fobj2,BestPop2,MeanPop2);
[BestCost3,CostB3,FEs3,Iter3,LayerNumber3,nPop3,Pop3,Cost3,FotonRate3,temp3,Costave3,UB3,LB3,VarNumber3,VarMin3,VarMax3,fobj3,BestPop3,MeanPop3]=CACS(FEs3,Iter3,LayerNumber3,nPop3,Pop3,Cost3,FotonRate3,temp3,Costave3,UB3,LB2,VarNumber3,VarMin3,VarMax3,fobj3,BestPop3,MeanPop3);
[BestCost4,CostB4,FEs4,Iter4,LayerNumber4,nPop4,Pop4,Cost4,FotonRate4,temp4,Costave4,UB4,LB4,VarNumber4,VarMin4,VarMax4,fobj4,BestPop4,MeanPop4]=CACS(FEs4,Iter4,LayerNumber4,nPop4,Pop4,Cost4,FotonRate4,temp4,Costave4,UB4,LB3,VarNumber4,VarMin4,VarMax4,fobj4,BestPop4,MeanPop4);
end
