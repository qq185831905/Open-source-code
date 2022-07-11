
function [BestCost,CostB,FEs,Iter,LayerNumber,nPop,Pop,Cost,FotonRate,temp,Costave,UB,LB,VarNumber,VarMin,VarMax,fobj,BestPop,MeanPop]=CACS(FEs,Iter,LayerNumber,nPop,Pop,Cost,FotonRate,temp,Costave,UB,LB,VarNumber,VarMin,VarMax,fobj,BestPop,MeanPop)
          
    Iter=Iter+1;    
    PopC=[];  
    CostC=[];
    NorDispCal=[];
    % Creat Quantum Layers
    MaxLay = randi(LayerNumber);
    NorDispInptut = 1:1:MaxLay;
    mu = 0;
    sigma = MaxLay/6;
    pd = makedist('Normal','mu',mu,'sigma',sigma); %创建一个正态分布对象，其均值mu = 0，标准差sigma = MaxLay/6。  
    NorDisp = pdf(pd,NorDispInptut);%y = pdf(gm,X)返回概率分布pd的概率密度函数(pdf)，在X的值处求值。  
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
                for ck=1:Dimension   %缩小
                    chaos(ck,1) = (Gbest_position(ck)-Fun_Lb)/(Fun_Ub-Fun_Lb); %只针对全局最优粒子进行混沌搜索,缩小成混沌变量
                end
                for ci=2:chaos_itrtn  %混沌迭代次数
                    for cj=1:Dimension   %遍历粒子的所有变量
                        chaos(cj,ci)= control * chaos(cj, ci - 1) * (1 - chaos(cj, ci - 1)); %logistic映射
                        tempchaos(cj,1)=(Fun_Ub-Fun_Lb) * chaos(cj,ci) + Fun_Lb;  %放大回原解空间
                    end
                    if chaos_constraint(tempchaos,Fun_Ub,Fun_Lb)==1   %是否符合限制条件
                        chaos_fitness(ci)=fobj(tempchaos');
                        if chaos_fitness(ci)<Gbest_Fitness  %替换
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

function y=chaos_constraint(tempchaos,Fun_Ub,Fun_Lb)   %混沌算法限制条件
a=max(tempchaos);
b=min(tempchaos);
if a<=Fun_Ub&&b>Fun_Lb
    y=1;
else
    y=0;
end
end