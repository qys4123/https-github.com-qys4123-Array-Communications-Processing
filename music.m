function SP = music(array, Rxx)

    derad = pi/180;        % deg -> rad
    radeg = 180/pi;
    h=size(array);
    N=h(1);                % N=5 the number of array
    M=h(2);                % M=3 the number of receiver
    d=min(array(:,1)):1:max(array(:,1));
    %d=0:1:4;
    
    [EV,D]=eig(Rxx);        % EV->eigenvectors, D->eigvalues
    EVA=diag(D)';           % compute the eigenvalue
    [EVA,I]=sort(EVA);      % from smallest to biggest 
    EV=fliplr(EV(:,I));     % from biggest to smallest(daoxu)
    En=EV(:,M+1:N);         % subspace of signal: 4-5column

for i=1:181
    
    angle=[i,0];
    Sp=spv(array, angle);
    Pmusic(i)=Sp'*En*En'*Sp;
    
end
    Pmusic=abs(Pmusic);
    SP=10*log10(Pmusic);
    
end