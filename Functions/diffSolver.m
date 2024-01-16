%diffSolver.m Solves the differential equations specified in Model.m either
%deterministically or stochastically

function [P, M, t_step, CT, MM,HILL]=diffSolver(solverParams)

P                 = solverParams.P;
M                 = solverParams.M;
Nt                = solverParams.Nt;
TauND_step        = solverParams.TauND_step; 
TauH_step         = solverParams.TauH_step;  
NM                = solverParams.NM;
gamma             = solverParams.gamma;
dm                = solverParams.dm;
dp                = solverParams.dp;      
dm1               = solverParams.dm1;  
dm2               = solverParams.dm2;
dp1               = solverParams.dp1; 
dp2               = solverParams.dp2;
dt                = solverParams.dt;
Stochastic        = solverParams.Stochastic;
rows              = solverParams.rows;
cols              = solverParams.cols;
percentCounter    = solverParams.percentCounter;
Up                = solverParams.Up;
Low               = solverParams.Low;
P_0Up             = solverParams.P_0Up;
n_Up              = solverParams.n_Up;
min_Up            = solverParams.min_Up;

cells=rows*cols;
CT=[(1:cells)' ,zeros(cells,Nt)];
MM=P.*0;
for t_step=1:Nt
    %% Percentage counter for the terminal
    if percentCounter==1
        if ~mod(t_step,floor(Nt*5/100)) 
            clc
            fprintf('Main loop progress: %.0f%% \n', t_step/Nt*100);
        end
    end
    

    
    %% Conditions for t<Tau (time delays)
    if max(TauND_step)+1>t_step
        Psum_delay=NM*P(:,1);
    else  
        if numel(TauND_step)==1
            Psum_delay=NM*P(:,t_step-TauND_step); %Average effect of Hes expressed by neighbouring 6 cells
        else
            Psum_delay=NM*P(sub2ind(size(P),[1:size(P,1)]',t_step-TauND_step));
        end
    end

    if max(TauH_step)+1>t_step
        p_delay=P(:,1);
    else
        if numel(TauH_step)==1
            p_delay=P(:,t_step-TauH_step);
        else
            p_delay=P(sub2ind(size(P),[1:size(P,1)]',t_step-TauH_step)); %Written this way to deal with the case when time delay is different across cells
        end
    end

    if min_Up==1
        a=Low;
        b=Up-a;
    else
        Up=(1-min_Up)./(1+(P_0Up./P(:,t_step)).^n_Up)+min_Up;
        a=Low;
        b=Up-a;
    end
    
    if Stochastic==1
      
        [new1, new2]=EulerStoch(cells,M(:,t_step),P(:,t_step),p_delay,Psum_delay,gamma,a,b,dm1,dm2,dp1,dp2,dt);
    else
        [new1, new2]=RK(M(:,t_step),P(:,t_step),p_delay,Psum_delay,gamma,a,b,dm,dp,dt);                  % 4th order Runge-Kutta solver
    end

    HILL(t_step)=1;

    M(:,t_step+1)=new1; % Storing new mRNA value
    P(:,t_step+1)=new2; % Storing new protein value


end
end


function [P_swap,P_binary]=probSwap(r,c,thresh)
    
    P_vec=rand(r*(c-1),1);
    P_swap=sparse(diag(P_vec,r));
    
    P_binary=P_swap;
    P_binary(P_swap>thresh)=1;
    P_binary(P_swap<thresh)=0;
    P_binary=sparse(P_binary);
    
%     figure(8667)
%     clf
%     h=surf(P_binary)
%     set(h, 'edgecolor','none')
%     view(0,90)
%     caxis([0 1])
%     colorbar
%     drawnow
end

function [new1, new2]=EulerStoch(cells,x1,x2,x3,x4,x5,x6,x7,f11,f12,f21,f22,dt)

    new1=x1 + f11(x1,x2,x3,x4,x5,x6,x7)*dt + f12(x1,x2,x3,x4,x5,x6,x7)*sqrt(dt).*normrnd(0,1,[cells, 1]);
    new2=x2 + f21(x1,x2,x3,x4,x5,x6,x7)*dt + f22(x1,x2,x3,x4,x5,x6,x7)*sqrt(dt).*normrnd(0,1,[cells, 1]);
%     new1=x1 + omega*f11(x1,x2,x3,x4,x5,x6,x7)*dt + sqrt(omega)*f12(x1,x2,x3,x4,x5,x6,x7)*sqrt(dt).*normrnd(0,1,[cells, 1]);
%     new2=x2 + omega*f21(x1,x2,x3,x4,x5,x6,x7)*dt + sqrt(omega)*f22(x1,x2,x3,x4,x5,x6,x7)*sqrt(dt).*normrnd(0,1,[cells, 1]);

    if min(new1)<0
        new1(new1<0)=0;
    end

    if min(new2)<0
        new2(new2<0)=0;
    end
end

%RK.m Implementation of the 4th-order Runge Kutta method. Outputs 
%the next step values for t+dt.

function [new1,new2]=RK(x1,x2,x3,x4,x5,x6,x7,f1,f2,dt)

    % x1 = Protein conc
    k1_1=f1(x1,       x2,x3,x4,x5,x6,x7)*dt;
    k2_1=f1(x1+k1_1/2,x2,x3,x4,x5,x6,x7)*dt;
    k3_1=f1(x1+k2_1/2,x2,x3,x4,x5,x6,x7)*dt;
    k4_1=f1(x1+k3_1,  x2,x3,x4,x5,x6,x7)*dt;

    % x2 = mRNA conc
    k1_2=f2(x1,x2,       x3,x4,x5,x6,x7)*dt;
    k2_2=f2(x1,x2+k1_2/2,x3,x4,x5,x6,x7)*dt;
    k3_2=f2(x1,x2+k2_2/2,x3,x4,x5,x6,x7)*dt;
    k4_2=f2(x1,x2+k3_2,  x3,x4,x5,x6,x7)*dt;

    %% New values
    new1=x1+(k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
    new2=x2+(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;

end

%swapRow.m Swaps rows I (vector of row numbers) with rows J (also vector of
%row numbers to be swapped)
function A=swapRow(A,I,J)
A([J,I],:)=A([I,J],:);
end