%% Modèle linéaire

%% Génération des données
NbRepetitionExperience=10;

for d=1:5;%1:5  %boucle sur la dimension
    
    for numMC=1:NbRepetitionExperience
    
    clear X;
    clear Y;
    clear Xhat;
    clear X_1;
    clear Xtilde;
    clear Wtilde;
    %d=2; %Dimension
    K=100; %Nombre de données dans l'échantillon
    for i=1:d
        X_1(i)=randn;
    end

    A_k=0.9.*eye(d);
    H_k=eye(d);

    sigmav=1;
    sigmaw=0.1;

    V=sigmav^2.*randn(d,K);
    W=sigmaw^2.*randn(d,K);

    X(:,1)=X_1;
    for k=2:K
        X(:,k)=A_k*X(:,k-1)+V(:,k);
        Y(:,k)=H_k*X(:,k)+W(:,k);
    end




    % Filtrage particulaire

    N=100; %Nombre de particules
    resample=0;
    Xtilde=zeros(K,N,d);
    Wtilde=1/N*ones(K,N,1); %initialisation des poids à 1/N
    Xtilde(1,:,:)=randn(1,N,d); %initialisation des particules

        for t=2:K
            for i=1:N
                Xtilde(t,i,:)=sigmav*randn(d,1)+A_k*reshape(Xtilde(t-1,i,:),1, d)';
                mat = (Y(:,t)-H_k*reshape(Xtilde(t,i,:),1,d)');
                Wtilde(t,i,1)=Wtilde(t-1,i,1)*exp((-(mat'*inv(sigmaw*eye(d))*mat)/2));
            end


    %         W=reshape(sum(Wtilde(t,:,:)),1,d);
    %         for i=1:N 
    %             Wtilde(t,i,:)=reshape(Wtilde(t,i,:),1,d)./W;
    %         end
            Wtilde(t,:,1)=Wtilde(t,:,1)/sum(Wtilde(t,:,1));


            Xhat(t,:)= sum(repmat(Wtilde(t,:,1)',1,d).*reshape(Xtilde(t,:,:),N,d));



            if (1/sum(Wtilde(t,:,1).^2))<N/2 %On peut jouer sur le paramètre N/2
                Xtilde_r=zeros(N,d);
             for i=1:N
                 u=rand;
                 j=1;
                 while sum(Wtilde(t,1:j,1))<u
                     j=j+1;
                 end

                 Xtilde_r(i,:)=Xtilde(t,j,:);
             end
             Xtilde(t,:,:)=Xtilde_r;    

            Wtilde(t,:)=1/N*ones(1,N,1); 
            end  


        end  

     Err(d,numMC)=sum(sum((X'-Xhat).^2));   
     figure(10);clf;plot(X');hold on;plot(Xhat,'--');%pause
     


    end
    X(:,1)=X_1;
for k=2:K
    X(:,k)=A_k*X(:,k-1)+V(:,k);
    Y(:,k)=H_k*X(:,k)+W(:,k);
end


end

%%
figure(13)
clf
plot(mean(Err,2));
title('Erreur quadratique en fonction de la dimension');
hold on
plot(mean(Err2,2), 'red');
   