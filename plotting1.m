% A Code to generate IMP rmse curve on Netflix Data Matrix 1

%% Part 1) Measure average RMSE of IMP
clear all

% Load Netflix Data Matrix 1
load data1
J = R>0; 
num = [1:5 10 20 30];

% Run IMP
for i = 1:length(num)
    
    n1 = num(i); n = n1*size(R,1);                                
    num_simulations = 20;

    gu = 8; gm = 8; 
    max_iter = 20;
     
    for sim_no = 1:num_simulations
        
        Je = J-Jt; 
        [a b] = find(Je==1); 
        tmp = sub2ind(size(Je),a,b); 
        k = randperm(length(tmp)); 
        indt = tmp(k(1:n));
        clear tmp k
        
        Je = zeros(size(Je)); Je(indt) = 1; 
        Re = Je.*R; Je = (Re>0); Rt = Jt.*R; 
                
        [Pr] = imp1(Re,Je,gu,gm);   
        [Pt Ps] = imp2(full(Re),Pr,gu,gm,max_iter);
        
        error(sim_no) = rmse(Rt,Jt,Pr,Ps,Pt,gu,gm);
        clear Re Je P* 
        
    end
    
    IMP(n1) = mean(error);
    
end

%% Part2) Plot RMSE curve
close all

figure
set(gca,'FontSize',18); 
plot(num,IMP(num),'s-r','LineWidth',3.5); 
axis([1,num(end),.9,1.40])
set(gca,'xtick',0:5:num(end));
set(gca,'ytick',.9:.05:1.4);
xlabel('Average Number of Observed Ratings per User')
ylabel('RMSE'); title('Netflix Data Matrix 1')
legend('IMP')
grid on
