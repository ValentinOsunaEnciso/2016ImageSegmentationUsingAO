%% Allostatic Optimization,  
%  Febrero, 2013, Valentin Osuna-Enciso, CUCEI-UDG, CIC-IPN
%  Evolve a Mixture of Gaussian functions 
%  with an automatic class number
function thres=AO_MGF
    clear all, format long
    %7Q:190,
    %6Q:55067,
    %5Q:302003,108073,198004,176035
    %4Q:388016,61060,120,296059,353013,253036,
    %3Q:250,159091,176039,118035,249-257,12074,66075,100080,106025,Im001_1
    %2Q:374067,35070,
    DB=imread('Im003_1.jpg');%
    n=3; D=n*3;         %Class number range
    DB=rgb2gray(DB);    %Convierto a escala de grises
    H=imhist(DB);       %Calculo del Histograma
    H=H/sum(H);         %Se normaliza Histograma experimental(suma de Hi=1)
    L=size(H,1);
    Amax=[max(H),L-1,L/(n*4)]; %Mmax=L-1;Dmax=L-1; 
    Amin=[0,1,0];
    x_h=[];x_l=[];
    for ind1=1:3
        x_h=[x_h;ones(n,1)*Amax(ind1)];
        x_l=[x_l;ones(n,1)*Amin(ind1)];
    end
    xp=0:1:255;
    % rng('shuffle'); %corridas=[];
    %%% 0 Initialization of parameters:
    Np=50; maxIter=1000; k=1; psi=0.03; %Q=randi(n-Qmin,Np,1);  Q=Q+Qmin-1; 
    %D=Q*3;%Dimensiones de cada individuo  
    %%% 1 Generate Home and Homeostasised Population:    
    Population=zeros(n*3,Np); %P2=Population;
    %% 2 Determine Allostatic Adaptation:
    for ind1=1:Np
        %Di=Q(ind1,1)*3;
        for ind2=1:D
            Population(ind2,ind1)=x_l(ind2,1)+rand()*(x_h(ind2,1)-x_l(ind2,1));  
        end
        Allo(ind1)=MGF(Population(:,ind1),H,n);
    end
    P2=Population;
    %Population=[Population;Q'];
    Best=Population(:,1); %D_best=Population(end,1);;%Best initial individual of population
    AlloBest=Allo(1); EvaluacioneFO=Np; cont=1; error(ind1)=999;%error,mov,temp
    for ind1=1:n*3
        HomeRange(ind1,1)=mean(Population(ind1,1:Np));
    end
    while AlloBest>0.12 && k<=maxIter   %EvaluacioneFO<=maxEvalFO% abs(0.00000001) 
      Stress1=randperm(Np,Np); otro=1;  
      for ind1=1:Np           
            %% 3 Create Mutant Individual,by only modifying 1 d:           
            SI=Population(:,ind1);
            Stress2=randi(D);
            SI(Stress2,1)=Population(Stress2,Stress1(1,ind1)); %Mutant Home Change         
            for ind2=1:n*3
                if(SI(ind2,1)<x_l(ind2,1))              
                    SI(ind2,1)=x_l(ind2,1);
                end
                if(SI(ind2,1)>x_h(ind2,1))          
                    SI(ind2,1)=x_h(ind2,1);
                end
            end
            %% 4 Determine Allostatic Adaptation between Individual&Mutant:                
            temp(1,ind1)=MGF(SI(:,1),H,n); 
            EvaluacioneFO=EvaluacioneFO+1;
            if temp(1,ind1)<Allo(1,ind1)
                 Allo(1,ind1)=temp(1,ind1); Population(:,ind1)=SI(:,1);  
                P2(:,ind1)=HomeRange; P2(Stress2,ind1)=SI(Stress2,1);               
                if temp(1,ind1)<AlloBest
                    otro=0; error(ind1)=AlloBest-temp(1,ind1);
                    AlloBest=temp(1,ind1); Best=SI;
                    mov(ind1)=(psi*(1-(psi/exp(psi*error(ind1)))))/(D);
                    mov2=mov(ind1)*2;                     
                    for ind3=1*round(Np/5)+1:round(Np/5)*2 
                       P2(:,ind3)=HomeRange(:,1)-mov(ind1)+mov2.*rand(D,1);
                       Population(:,ind3)=HomeRange(:,1)+HomeRange(:,1).*randn(D,1)*mov(ind1);
                    end
                    for ind3=2*round(Np/5)+1:round(Np/5)*3
                        P2(:,ind3)=Best(:,1)-mov(ind1)+mov2.*rand(D,1);
                        Population(:,ind3)=Best(:,1)+Best(:,1).*randn(D,1)*mov(ind1);
                    end
                    r1=max(Population(1:end,:)');
                    r2=min(Population(1:end,:)');r2=r2';r1=r1';
                    for ind3=3*round(Np/5)+1:Np %Sigue explorando   
                        P2(:,ind3)=(HomeRange(:,1).*rand(D,1));
                        Population(:,ind3)=r2+(r1-r2).*rand(D,1);
                    end
                    Population(Stress2,Stress1(1,ind1))=SI(Stress2,1);
                end
            else
                %% 5 Modify Allostatic Load if there weren't adaptation:
%                 r1=max(P2(1:end,:)');r2=min(P2(1:end,:)');r2=r2';r1=r1';
%                     for ind3=3*round(Np/5)+1:Np %Sigue explorando   
%                         P2(:,ind3)=(HomeRange(:,1).*rand(D,1));
%                         Population(:,ind3)=r2+(r1-r2).*rand(D,1);
%                     end
            end               
      end %for ind1=1:Np
      %% 6 Modifying population:
      if otro==1
        for ind1=1:D
            HomeRange(ind1,1)=mean(P2(ind1,1:Np));
        end
      end
      fprintf('k=%d,f(best)=%.6f\n',k,AlloBest);
      CONVERGENCIA(cont)=AlloBest; cont=cont+1; 
      k=k+1; 
  end %while AllostaticAdaptation(1,1)>limit && k<maxIter 
  thres=grafica(Best,H,DB,D,n);
end %function