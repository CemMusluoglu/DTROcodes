function [dist,path]=shortest_path(adj,q)
    K=size(adj,1);
    dist=Inf(K,1);
    dist(q)=0;
   
    visited=[];
    pred=zeros(1,K);
    unvisited=setdiff([1:K],visited);
    path=cell(1,K);

    while(length(visited)<K)
        I=find(dist==min(dist(unvisited)));
        I=I';
       
        for ind=I
            visited=[visited,ind];
            unvisited=setdiff([1:K],visited);
            neighbors_i=find(adj(ind,:)==1);
            for m=intersect(neighbors_i,unvisited)
                if(dist(ind)+1<dist(m))
                    dist(m)=dist(ind)+1;
                    pred(m)=ind;
                end
            end
        end
        
    end
    
    for k=1:K
        jmp=k;
        path_k=[k];
        while(jmp~=1)
            jmp=pred(jmp);
            path_k=[jmp,path_k];
        end
        path(k)={path_k};
    end
    
end