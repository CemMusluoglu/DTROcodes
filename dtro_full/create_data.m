function [Ufull,Vfull]=create_data(nbsensnode,nbsamples,nbnodes,Q)


D=10;    %dimension of latent random process
noisepower=0.1; %additive spatially uncorrelated white noise with uniform distribution
signalvar=0.5;

rng('shuffle');
d=randn(nbsamples,Q); %random latent process
d=sqrt(signalvar)./(sqrt(var(d))).*(d-ones(nbsamples,1)*mean(d));
s=randn(nbsamples,D-Q); %same random latent process, but new observations
s=sqrt(signalvar)./(sqrt(var(s))).*(s-ones(nbsamples,1)*mean(s));

for k=1:nbnodes
    Ainit{k}=rand(Q,nbsensnode(k))-0.5; %random mixture
    Binit{k}=rand(D-Q,nbsensnode(k))-0.5; %random mixture
    noise{k}=(randn(nbsamples,nbsensnode(k))); %same random noise, but new observations
    noise{k}=sqrt(noisepower)./sqrt(var(noise{k})).*(noise{k}-ones(nbsamples,1)*mean(noise{k}));
end

teller=0;

for k=1:nbnodes
    V{k}=s*Binit{k}+noise{k};
    U{k}=d*Ainit{k}+V{k};

    Ufull(1:nbsamples,teller+1:teller+nbsensnode(k))=U{k};
    Vfull(1:nbsamples,teller+1:teller+nbsensnode(k))=V{k};
    teller=teller+nbsensnode(k);
end

end