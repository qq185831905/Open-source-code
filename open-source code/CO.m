function [Gbest_position,Gbest_Fitness]=CO(Fun_Ub,Fun_Lb,Gbest_position,Gbest_Fitness,fobj)
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
           chaos_fitness=fobj(tempchaos);
        if chaos_fitness(ci)<Gbest_Fitness  %替换
            Gbest_Fitness=chaos_fitness(ci);
            Gbest_position=tempchaos;
      
        end
    end
end
end

%初始化：
function y=chaos_constraint(tempchaos,Fun_Ub,Fun_Lb)   %混沌算法限制条件
a=max(tempchaos);
b=min(tempchaos);
if a<=Fun_Ub&&b>Fun_Lb
    y=1;
else
    y=0;
end
end

