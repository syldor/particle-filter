% Author: Sylvain Dorey
% Contact: sylvain.dorey@centraliens-lille.org
%
% Some simulations of the particle filter

%% Filtrage particulaire
clear Xhat
clear X
N=100; %Number of particules
T=0.01;
data=load('data.mat');
X=data.x;
Y=data.y;
resample=0;
Xtilde=zeros(501,N);
Wtilde=1/N*ones(501,N); % Initialization of weights at 1/N
Xtilde(1,:)=sqrt(10)*randn(1,N); % Initialization of particules
%%
for t=2:501
    for i=1:N
        Xtilde(t,i)=sqrt(10)*randn+f(Xtilde(t-1,i),(t-1)*T);
        Wtilde(t,i)=Wtilde(t-1,i)*(1/sqrt(2*pi)*exp((-(Y(t)-g(Xtilde(t,i)))^2)/2));
    end
    W=sum(Wtilde(t,:));
    for i=1:N 
        Wtilde(t,i)=Wtilde(t,i)/W;
    end
    
    Xhat(t)= Wtilde(t,:)*Xtilde(t,:)';  
     
    % Multinomial resampling
    %figure(2);
    %plot(Wtilde(t,:))
    %title(num2str([1/sum(Wtilde(t,:).^2)]))
    %pause
    
    if (1/sum(Wtilde(t,:).^2))<N/2 % We can play with the parameter N/2
        resample=resample+1;     % To count the number of resamplings
     for i=1:N
         u=rand;
         j=1;
         while sum(Wtilde(t,1:j))<u
             j=j+1;
         end
         
         Xtilde_r(i)=Xtilde(t,j);
     end
     Xtilde(t,:)=Xtilde_r;    
  
    Wtilde(t,:)=1/N*ones(1,N); 
    end  
    
    

end

Err=sum((X-Xhat).^2) 
%%
% Display of results
figure(1)
clf
plot(X(1:size(Xhat,2)))
hold on
plot(Xhat,'r')
title('Filtre particulaire')
legend('Valeur exacte','Valeur estimée')

%%
%% Graph of the weights
% We look at each weight one after one, in a graph (particle, weight),
% press enter to switch.


figure(2)
for t=1:100;
%clf
plot(Wtilde(t,:))
hold on;
%pause();
end
%%
figure(123)
AxisPDF=[-10:0.1:10];
for t=1:6;
subplot(2,3,t);
ApproxPosteriorPDF=ksdensity(Xtilde(t,:),AxisPDF,'weights',Wtilde(t,:)); 
figure(10);clf;plot(AxisPDF,ApproxPosteriorPDF);
hold on,;plot(X(t),max(ApproxPosteriorPDF),'or')
hold on,;plot(Y(t),max(ApproxPosteriorPDF),'xg')
hold on,;plot(Xhat(t),max(ApproxPosteriorPDF),'sk')
legend('Distribution des poids', 'Valeur réelle', 'Valeur observée', 'Valeur estimée');
pause
end




%%
figure(3)
plot(Xtilde(5,:),Wtilde(5,:),'.')

%% Linear model

%% Data generation

%for d=1:20  %loop on the dimension
clear X;
clear Y;
clear Xhat;
clear X_1;
d=1; % Dimension
K=100; % Number of data in the sample
for i=1:d
    X_1(i)=randn;
end

A_k=0.9.*eye(d);
H_k=eye(d);

sigmav=1;
sigmaw=1;

V=sigmav^2.*randn(d,K);
W=sigmaw^2.*randn(d,K);

X(:,1)=X_1;
for k=2:K
    X(:,k)=A_k*X(:,k-1)+V(:,k);
    Y(:,k)=H_k*X(:,k)+W(:,k);
end


%%

% Filtrage particulaire

N=100; 
resample=0;
Xtilde=zeros(K,N,d);
Wtilde=1/N*ones(K,N,d); 
Xtilde(1,:,:)=randn(1,N,d); 

    for t=2:K
        for i=1:1
            Xtilde(t,i,:)=sigmav*randn(d,1)+A_k*reshape(Xtilde(t-1,i,:),1, d)';
            mat = (Y(:,t)-H_k*reshape(Xtilde(t,i,:),1,d)');
            Wtilde(t,i,:)=reshape(Wtilde(t-1,i,:),1,d)*inv((sigmaw*eye(d)*sqrt(2*pi)))*exp((-(mat'*inv(sigmaw*eye(d))*mat)/2));
        end


        W=reshape(sum(Wtilde(t,:,:)),1,d);
        for i=1:N 
            Wtilde(t,i,:)=reshape(Wtilde(t,i,:),1,d)./W;
        end


        Xhat(t,:)= reshape(Wtilde(t,:,1),100,1)'*reshape(Xtilde(t,:,:),100,d);
    end  
    if (1/sum(Wtilde(t,:,1).^2))<N/2 
     for i=1:N
         u=rand;
         j=1;
         while sum(Wtilde(t,1:j,1))<u
             j=j+1;
         end
         
         Xtilde_r(i)=Xtilde(t,j,:);
     end
     Xtilde(t,:,:)=Xtilde_r;    
  
    Wtilde(t,:)=1/N*ones(1,N,d); 
    end  

    
 Err(d)=sum(sum((X'-Xhat).^2));
end 
%%
figure(13)
clf
plot(Err);    
    

%%

 


figure(15)
clf
plot(X(1:size(Xhat,2)))
hold on
plot(Xhat,'r')
hold on
plot(Y(1:size(Xhat,2)), 'black')
title('Filtre particulaire')
legend('Valeur exacte','Valeur estimée')

%%
%% Génération des données
for num = 1:10;
for d=1:5  %boucle sur la dimension
clear X;
clear Y;
clear Xhat;
clear X_1;
%d=1; %Dimension
K=100; %Nombre de données dans l'échantillon
for i=1:d
    X_1(i)=randn;
end

A_k=0.9.*eye(d);
H_k=eye(d);

sigmav=1;
sigmaw=1;

V=sigmav^2.*randn(d,K);
W=sigmaw^2.*randn(d,K);

X(:,1)=X_1;
for k=2:K
    X(:,k)=A_k*X(:,k-1)+V(:,k);
    Y(:,k)=H_k*X(:,k)+W(:,k);
end

%Kalman
%initialisation

Xhat(:,1) = Y(:,1);
Pk1k1 = eye(d);
for k = 2:K
    %Prediction
    Xkk1=A_k*Xhat(:,k-1);
    Pkk1=A_k*Pk1k1*A_k'+sigmav*eye(d);

    %Mise à jour
    ytildek=Y(:,k)-H_k*Xkk1;
    Sk=H_k*Pkk1*H_k'+sigmaw*eye(d);
    Kk=Pkk1*H_k'*inv(Sk);
    Xhat(:,k)=Xkk1+Kk*ytildek;
    Pkk=(eye(d) - Kk*H_k)*Pkk1;

end
Err2(d,num)=sum(sum((X-Xhat).^2));
end
end

%%
figure(145)
plot(X);
hold on :
plot(Xhat, 'r');

%%
figure(33)
plot(mean(Err2,2))

