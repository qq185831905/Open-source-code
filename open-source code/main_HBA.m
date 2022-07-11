%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Honey Badger Algorithm source code 
%     "Honey Badger Algorithm: New Metaheuristic Algorithm for %  %     Solving Optimization Problems." 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;
close all;
fitfun = @sumsqu;
dim=30; %维度
T=1000; %迭代次数
Lb=-10; %下限
Ub=10; %上限
N=30;
[xmin,fmin,CNVG]=HBA(fitfun,dim,Lb,Ub,T,N);
figure,
semilogy(CNVG,'r')
xlim([0 T]);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
legend('HBA')

% display(['The best location= ', num2str(xmin)]);
% display(['The best fitness score = ', num2str(fmin)]);

        

