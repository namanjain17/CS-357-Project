
pkg load statistics;

function i=RouletteWheelSelection(P)

    r=rand;
    C=cumsum(P);
    i=find(r<=C,1,'first');

end

nOfConsEq = 0;
nOfConsIneq = 12;
g = @(x) [-x(1)+0.0193*x(3),-x(2)+0.00954,-pi*x(3)^2*x(4)-pi*x(3)^2/3+1296000,x(4)-240,-x(1),-x(2),x(1)-99,x(2)-99,10-x(3),10-x(4),x(3)-200,x(4)-200];
h = @(x) [0];

MR=0.9;

CostFunction=@(x) (0.6244*x(1)*x(3)*x(4) + 1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3));        % Cost Function


function i=CalculateDebb(x,y)
        CostFunction=@(x) (0.6244*x(1)*x(3)*x(4) + 1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3));
        g = @(x) [-x(1)+0.0193*x(3),-x(2)+0.00954,-pi*x(3)^2*x(4)-pi*x(3)^2/3+1296000,x(4)-240,-x(1),-x(2),x(1)-99,x(2)-99,10-x(3),10-x(4),x(3)-200,x(4)-200];
        h = @(x) [0];
        temp_costx = CostFunction(x);
        temp_costy = CostFunction(y);

        Gx=g(x);
        Hx=h(x);
        Gy=g(y);
        Hy=h(y);

        Ax = Gx>0;
        Bx = Hx>0;
        temp_feasiblex = ((sum(Ax)+sum(Bx))==0);
        Ay = Gy>0;
        By = Hy>0;
        temp_feasibley = ((sum(Ay)+sum(By))==0);
        
        if(temp_feasiblex==0 && temp_feasibley==1) 
             i=0;
        elseif(temp_feasiblex==1 && temp_feasibley==0)
              i=1;
        elseif(temp_feasiblex==1 && temp_feasibley==1)
              i=(temp_costx<=temp_costy); 
        else
             i=(sum(Gx(Ax))+sum(Hx(Bx)))<=(sum(Gy(Ay))+sum(Hy(By)));
        end
end




nVar=4;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=[0,0,10,10];         % Decision Variables Lower Bound
VarMax=[99,99,200,200];         % Decision Variables Upper Bound




%% ABC Settings

MCN=1000;              % Maximum Number of Iterations

SN=250;               % Population Size (Colony Size)

nOnlooker=SN;         % Number of Onlooker Bees

Limit=20; % Abandonment Limit Parameter (Trial Limit)

a=1;                    % Acceleration Coefficient Upper Bound

One = ones(1,nVar);

%% Initialization

% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];
empty_bee.Debb =[];

% Initialize Population Array
pop=repmat(empty_bee,SN,1);

% Initialize Best Solution Ever Found
BestSol.Cost=Inf;
BestSol.Debb=-Inf;
BestSol.Position=[-Inf,-Inf,-Inf,-Inf];
% Create Initial Population


for i=1:SN
    for j=1:nVar
        pop(i).Position(j)=unifrnd(VarMin(j),VarMax(j),1);
    end

    pop(i).Cost=CostFunction(pop(i).Position);
    if (CalculateDebb(pop(i).Position,BestSol.Position))
        BestSol=pop(i);
    end
end




% Abandonment Counter
C=zeros(SN,1);

% Array to Hold Best Cost Values
BestCost=zeros(MCN,1);

%% ABC Main Loop


for it=1:MCN
    
    % Recruited Bees
    for i=1:SN
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:SN];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a.*unifrnd(0,+1,1,nVar);
        %j = randi([1:nVar]);
        R = unifrnd(0,1,1,nVar);
        R = R<MR;
        % New Bee Position
        newbee=pop(i);
        newbee.Position=pop(i).Position+(R.*phi).*(pop(i).Position-pop(k).Position);
        % if(it>=0.75*MCN) 
            newbee.Position-=(R.*(One-phi)).*(pop(i).Position-BestSol.Position);
       % end
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        newbee.Debb=0;
        % Comparision
        if (CalculateDebb(newbee.Position,pop(i).Position))
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end

    end


    % Calculate Fitness Values and Selection Probabilities

    
    F=zeros(SN,1);
    MeanCost = mean([pop.Cost]);
    feasible=zeros(1,SN);
    CV=zeros(1,SN);
   
    for i=1:SN
        F(i) = exp(-pop(i).Cost/ MeanCost);
       % if(pop(i).Cost>=0)
       %   F(i) = 1/(1+pop(i).Cost);
       %  else
       %   F(i) = 1+abs(pop(i).Cost);
       % end
        G=g(pop(i).Position);
        H=h(pop(i).Position);
        A = G>0;
        B = H>0;
        feasible(i) = ((sum(A)+sum(B))==0);
        CV(i) = (sum(G(A))+sum(H(B)));
    end
 
    Sum_F = sum(F);
    Sum_CV = sum(CV);

    P=zeros(SN,1);
    for i=1:SN
          if(feasible(i) == 1)
              P(i) = 0.5 + ((0.5*F(i))/Sum_F);
          else
              P(i) = (1-(CV(i)/Sum_CV))*0.5;
          end
    end
     
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:SN];
        k=K(randi([1 numel(K)]));
        
        % New Bee Position
        
        phi=a.*unifrnd(0,+1,1,nVar);
        %j = randi([1:nVar]);
        R = unifrnd(0,1,1,nVar);
        R = R<MR;
        % New Bee Position
        newbee=pop(i);
        newbee.Position=pop(i).Position+(R.*phi).*(pop(i).Position-pop(k).Position);
        % if(it>=0.75*MCN) 
            newbee.Position-=(R.*(One-phi)).*(pop(i).Position-BestSol.Position);
       % end

        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);


        % Comparision
        if (CalculateDebb(newbee.Position,pop(i).Position))
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:SN
        if C(i)>=Limit
          if(0)
             for j=1:nVar
                pop(i).Position(j)=unifrnd(VarMin(j),VarMax(j),1);
             end
           else
                phi=a.*unifrnd(0,+1,1,nVar);
                %xlim=unifrnd(VarMin(j),VarMax(j),1);
                pop(i).Position = pop(i).Position - phi.*(pop(i).Position-BestSol.Position);
           end
              pop(i).Cost=CostFunction(pop(i).Position);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:SN
        if (CalculateDebb(pop(i).Position,BestSol.Position))
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
    
%% Results

figure;
plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

