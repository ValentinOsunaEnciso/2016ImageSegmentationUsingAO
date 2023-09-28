%% FUNCION GRAFICA GAUSSIANAS E IMAGEN SEGMENTADA: %%%%%%%%%%%%%%%%%%%%%%%%
% Recibe x_best: La mejor particula, D: dimensiones de cada particula
% H: histograma de la imagen, DB: imagen en escala de gris
function T1=grafica(x_best,H,DB,D_best,n)
    xp=0:1:255; alturas=x_best(1:n,1)'; medias=x_best(n+1:(2*n),1)'; desvest=x_best(2*n+1:(3*n),1)';
%     for ind1=1:3:D_best
%        alturas=[alturas,x_best(ind1,1)]; 
%        medias=[medias,x_best(ind1+1,1)];
%        desvest=[desvest,x_best(ind1+2,1)];
%     end
    [medias ind]=sort(medias); temp=alturas(ind(1:end)); alturas=temp; 
    temp=desvest(ind(1:end)); desvest=temp; tam=size(medias,2);
    xp=0:1:255;mix=zeros(256,1);xp=xp'; %figure; hold on;
    for ind1=1:tam
            %%%Gauss
        mix=mix+(alturas(ind1)*exp(-((xp-round(medias(ind1))).^2)/...
            (2*(desvest(ind1)^2))));%Gauss
%         plot(alturas(ind1)*exp(-((xp-round(medias(ind1))).^2)/...
%           (2*(desvest(ind1)^2))),'g--','LineWidth',2);%Gauss
            %%%Laplacian
%         mix=mix+(alturas(ind1)*exp(-((xp-round(medias(ind1))))/...
%             (2*(desvest(ind1)^2))));%Laplacian
%         plot(alturas(ind1)*exp(-((xp-round(medias(ind1))))/...
%           ((desvest(ind1)))),'g--','LineWidth',2);%Gauss
            %%%Cauchi
% %         mix=mix+(alturas(ind1).*((desvest(ind1)^2)./(((xp-...
% %             round(medias(ind1))).^2)+(desvest(ind1)^2))));%Cauchi3parameter
%         plot(alturas(ind1).*((desvest(ind1)^2)./(((xp-...
%             round(medias(ind1))).^2)+(desvest(ind1)^2))),'k--','LineWidth',2);%Cauchi3parameter
            %%%Cauchy
%         mix=mix+alturas(ind1)./(1+(((xp-round(medias(ind1))).^2)./desvest(ind1)));
%         plot(alturas(ind1)./(1+(((xp-round(medias(ind1))).^2)./desvest(ind1))),'k--','LineWidth',2);
    end
%     plot(H,'r','LineWidth',2); figure; hold on;
%     plot(mix,'b','LineWidth',2);
%     plot(H,'r','LineWidth',2); 
    %%%%%       CALCULATE THRESHOLDS:                  %%%%%%%%%%%%%%%%%%%%
    T=[];Q=size(medias,2);
    for ind1=1:Q-1
        a=0; b=0; c=0;
        a=(desvest(1,ind1)^2)-(desvest(1,ind1+1)^2);        
        if a~=0            
            b=2*((medias(1,ind1)*(desvest(1,ind1+1)^2))-...
            (medias(1,ind1+1)*(desvest(1,ind1)^2)));
            c=((desvest(1,ind1)*medias(1,ind1+1))^2)-...
            ((desvest(1,ind1+1)*medias(1,ind1))^2)+...
            (2*((desvest(1,ind1)*desvest(1,ind1+1))^2)*...
            log((desvest(1,ind1+1)*alturas(1,ind1))/(desvest(1,ind1)*...
            alturas(1,ind1+1))));
            T1a=(-b+sqrt((b^2)-(4*a*c)))/(2*a);
            T1b=(-b-sqrt((b^2)-(4*a*c)))/(2*a);
            if(T1a<=medias(1,ind1+1) && T1a>=medias(1,ind1))  
                T=[T,T1a];
            else
                T=[T,T1b];
            end  
        else
            T1a=(medias(1,ind1+1)+medias(1,ind1))/2;
            T=[T,T1a];
        end
              
    end
    T1=sort(T);
    fprintf('T=%d \n',T);
    T=[0,T,255];
    %%%%%       PAINTING THE PIXELS:    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colorxgrupo=1/Q; inicial=0;
    [fila columna]=size(DB);
    for ind0=1:Q
      for ind1=1:fila
        for ind2=1:columna
            if (DB(ind1,ind2)<T(1,ind0+1))&&(DB(ind1,ind2)>=T(1,ind0))
                DBsegmented(ind1,ind2)=inicial;
            end
        end
      end
      inicial=inicial+colorxgrupo;
    end
    figure,DBsegmented=mat2gray(DBsegmented);
    imshow(DBsegmented)
    %% Modified Hausdorff distance:
    % direccion=['C:\Users\sossa2012x\Desktop\2013_AutomaticThreholding\Img\','Im238_0_GT','.tif'];
    % AI=DBsegmented>0;
    % BI=imread(direccion);BI=rgb2gray(BI(:,:,1:3));
    % [AIbordes t]=edge(AI,'canny',.1);
    % [BIbordes t]=edge(BI,'canny',.1);
    % [A(:,1) A(:,2)]=find(AIbordes);
    % [B(:,1) B(:,2)]=find(BIbordes);
    % [ mhd ] = ModHausdorffDist( A, B )
end