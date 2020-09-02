function [Ufull,Vfull]=create_data(nbsensnode,nbsamples,nbnodes,Q)


D=10;    %dimension of latent random process
noisepower=0.1; %additive spatially uncorrelated white noise with uniform distribution

rng('shuffle');
d1=randn(nbsamples,Q); %random latent process
d1=sqrt(1)*(d1-ones(nbsamples,1)*mean(d1));
d2=randn(nbsamples,D-Q); %same random latent process, but new observations
d2=sqrt(1)*(d2-ones(nbsamples,1)*mean(d2));

for k=1:nbnodes
    Ainit{k}=rand(Q,nbsensnode(k))-0.5; %random mixture
    Binit{k}=rand(D-Q,nbsensnode(k))-0.5; %random mixture
    noise1{k}=sqrt(noisepower)*(rand(nbsamples,nbsensnode(k))); %random noise
    noise1{k}=noise1{k}-ones(nbsamples,1)*mean(noise1{k});
    noise2{k}=sqrt(noisepower)*(randn(nbsamples,nbsensnode(k))); %same random noise, but new observations
    noise2{k}=noise2{k}-ones(nbsamples,1)*mean(noise2{k});
end

teller=0;

for k=1:nbnodes
    V{k}=d2*Binit{k}+noise2{k};
    U{k}=d1*Ainit{k}+V{k};

    Ufull(1:nbsamples,teller+1:teller+nbsensnode(k))=U{k};
    Vfull(1:nbsamples,teller+1:teller+nbsensnode(k))=V{k};
    teller=teller+nbsensnode(k);
end

end