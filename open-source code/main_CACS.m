
%%************ core code of CACS *************
%clear all
%mex cec14_func.cpp -DWINDOWS
clc
clear
func_num=1;

%�ܲ���
pop_size=10;
iter_max=500;

% Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Function_name='F23'; 
% Load details of the selected benchmark function
[Xmin,Xmax,D,fobj]=get_functions_details30(Function_name);
Max_FES=10000*D;
        
%�趨һ����ʼ����Ⱥ
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
%setpop ��ʼ��Ⱥ
%MaxFes Ĭ��10000*����ά��=������Ӧ�ȵĴ���
%VarNumber ����ά��
%nPop ��Ⱥ����
%T �������ֵ
%LB ����
%UB ����
%fobj Ŀ�꺯��
%% Problem Information
VarMin = LB *ones(1,VarNumber);        % Lower bound of variable;
VarMax = UB *ones(1,VarNumber);         % Upper bound of variable;
BestCost=10000;
%% General Parameters of the Algorithm

LayerNumber = 5 ;     % Maximum number of Layers around nucleus
FotonRate = 0.1 ;     % Initial Foton Rate for position determination of electrons

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
Pop=Pop(SortOrder,:); %��˳���ź�
BestPop=Pop(1,:);  
MeanPop=mean(Pop); %ÿ�еľ�ֵ

%% Main Loop
while FEs<MaxFes           
    Iter=Iter+1;    
    PopC=[];  
    CostC=[];
    NorDispCal=[];
    % Creat Quantum Layers
    MaxLay = randi(LayerNumber);
    NorDispInptut = 1:1:MaxLay;
    mu = 0;
    sigma = MaxLay/6;
    pd = makedist('Normal','mu',mu,'sigma',sigma); %����һ����̬�ֲ��������ֵmu = 0����׼��sigma = MaxLay/6��  
    NorDisp = pdf(pd,NorDispInptut);%y = pdf(gm,X)���ظ��ʷֲ�pd�ĸ����ܶȺ���(pdf)����X��ֵ����ֵ��  
    NorDispCal(1,:)=NorDisp;
    NorDispCal(2,:)=NorDispCal(1,:)./sum(NorDispCal(1,:));
    NorDispCal(3,:)=nPop*NorDispCal(2,:);
    NorDispCal(4,:)=round(NorDispCal(3,:));
    NorDispCal(5,:)=cumsum(NorDispCal(4,:));   
    LayCol=[0 NorDispCal(5,:)];
    LayCol(LayCol>nPop)=nPop;   
    % Search Loop   
    for i=1:MaxLay        
        PopA=Pop(LayCol(i)+1:LayCol(i+1) , :);
        CostA=Cost(LayCol(i)+1:LayCol(i+1) , 1);
        Energy=mean(CostA);
        Orbit=i;     
        for j=1:size(PopA,1)  
               
         
           if rand>FotonRate                
                if CostA(j,1)>Energy                    
                    Ir=unifrnd(0,1,1,2);
                    Jr=unifrnd(0,1,1,VarNumber);                   
                    Xold=PopA(j,:);
                    Xbest=BestPop;
                    Xmean=MeanPop;                   
                    PopB(j,:)=Xold+(Jr.*(Ir(1)*Xbest-Ir(2)*Xmean)/Orbit);                                         
                    PopB(j,:) = max(PopB(j,:),VarMin);
                    PopB(j,:) = min(PopB(j,:),VarMax);                    
                    CostB(j,1)=fobj(PopB(j,:));  
                       % 	Dynamic Photons
           if CostB(j,1)>Costave
               FotonRate=0.1+rand*0.1;
           else
               FotonRate=0.2-0.1*i/MaxLay;
           end
                    FEs=FEs+1;                   
                else                    
                    Ir=unifrnd(0,1,1,2);
                    Jr=unifrnd(0,1,1,VarNumber);                    
                    Xold=PopA(j,:);
                    Xbest=PopA(1,:);
                    if size(PopA,1)==1
                        Xmean=PopA; 
                    else
                        Xmean=mean(PopA);
                    end
                    PopB(j,:)=Xold+(Jr.*(Ir(1)*Xbest-Ir(2)*Xmean));
                    PopB(j,:) = max(PopB(j,:),VarMin);
                    PopB(j,:) = min(PopB(j,:),VarMax);                    
                    CostB(j,1)= fobj(PopB(j,:));                    
                    FEs=FEs+1;                    
                end                
            else     %Chaotic orbit   
                Fun_Ub=UB;
                Fun_Lb=LB;
                Gbest_position=PopA(j,:)';
                Gbest_Fitness=CostA(j,1);

                chaos_itrtn = 10;
                Dimension = size(Gbest_position,1);
                control=4;
                for ck=1:Dimension   %��С
                    chaos(ck,1) = (Gbest_position(ck)-Fun_Lb)/(Fun_Ub-Fun_Lb); %ֻ���ȫ���������ӽ��л�������,��С�ɻ������
                end
                for ci=2:chaos_itrtn  %�����������
                    for cj=1:Dimension   %�������ӵ����б���
                        chaos(cj,ci)= control * chaos(cj, ci - 1) * (1 - chaos(cj, ci - 1)); %logisticӳ��
                        tempchaos(cj,1)=(Fun_Ub-Fun_Lb) * chaos(cj,ci) + Fun_Lb;  %�Ŵ��ԭ��ռ�
                    end
                    if chaos_constraint(tempchaos,Fun_Ub,Fun_Lb)==1   %�Ƿ������������
                        chaos_fitness(ci)=fobj(tempchaos');
                        if chaos_fitness(ci)<Gbest_Fitness  %�滻
                            Gbest_Fitness=chaos_fitness(ci);
                            Gbest_position=tempchaos';

                        end
                    end
                end
                PopB(j,:)=Gbest_position;
                CostB(j,1)=Gbest_Fitness;
                temp=temp+Gbest_Fitness;
                Costave=temp/i;
                FEs=FEs+1;  
           end  

        end    
        PopC=[PopC ; PopB];
        CostC=[CostC ; CostB];        
    end
    
    % Merge Candidates
    Pop=[Pop ; PopC];
    Cost=[Cost ; CostC];   

    % Sort Population
    [Cost, SortOrder]=sort(Cost);
    Pop=Pop(SortOrder,:);
    BestPop=Pop(1,:);
    BestCost=Cost(1,1);
    MeanPop=mean(Pop);    
    Pop=Pop(1:nPop,:);
    Cost=Cost(1:nPop,:);
    
    % Store Best Cost Ever Found
    BestCosts(Iter,1)=BestCost;

    % Show Iteration Information
    disp(['Iteration ' num2str(Iter) ': Best Cost = ' num2str(BestCosts(Iter))]);
end

function y=chaos_constraint(tempchaos,Fun_Ub,Fun_Lb)   %�����㷨��������
a=max(tempchaos);
b=min(tempchaos);
if a<=Fun_Ub&&b>Fun_Lb
    y=1;
else
    y=0;
end
end





