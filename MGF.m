%% Mixture of Gaussian Functions: it works with 3 Gaussian functions. %%%%
function error=MGF(x,H,n)

xp=0:1:255;mix=zeros(256,1);xp=xp';
    for ind1=1:n
    mix=mix+(x(ind1,1)*exp(-((xp-round(x(ind1+n,1))).^2)/(2*(x(ind1+(n*2),1)^2))));%Gauss
%     mix=mix+x(ind1,1)./(1+(((xp-round(x(ind1+1,1))).^2)./x(ind1+2,1)));%Cauchy??
%     mix=mix+(x(ind1,1).*((x(ind1+(n*2),1)^2)./(((xp-round(x(ind1+n,1))).^2)+(x(ind1+(n*2),1)^2))));%Cauchi
%     mix=mix+(x(ind1,1)*exp(-((xp-round(x(ind1+1,1))))/((x(ind1+2,1)))));%Laplacian
    end    
    medias=[]; cont=0;
    for ind=n+1:n*2
        medias=[medias,x(ind,1)];
        cont=cont+1;
    end
    medias=sort(medias); suma=0;
    for ind=1:cont-1
        if(medias(1,ind+1)-medias(1,ind))<(256*1.5)/(n*3)
            suma=suma+(1/(medias(1,ind+1)-medias(1,ind)));
        end
    end
%     medias=[round(x(4,1)),round(x(5,1)),round(x(4,1))];
%     medias=sort(medias);
%     error=HellingerDistance(H,mix)+(3/((medias(3)-medias(2))+(medias(2)-medias(1))+eps));
    %error=HellingerDistance(H,mix)+(suma/cont)+(Di/100);
    error=HellingerDistance(H,mix)+(suma);
    %error=HellingerDistance(H,mix)+(Di/150);
%     error=HellingerDistance(H,mix);
end
%% 
    %[M temp]=sort(round(x(4:6,1)));
    %penalizacion=w*abs(sum(mix)-1);            % Original
    %penalizacion=1/((M(3)-M(1))+(M(2)-M(1)));  % Modificada
    % Distancia Euclideana:
%     error=((H-mix)).^2;
%     error=sum(error)/size(H,1);
%     error=error+penalizacion;
%     error1=MahalanobisDistance(H,mix);
%     error2=BhattacharyyaDistance(H,mix);
    % Distancia Minkowski:
%   error4=abs(H-mix);[error4 temp]=max(error4);error4=error4+penalizacion;
%     error4=abs(H-mix);error4=sum(error4)+penalizacion;