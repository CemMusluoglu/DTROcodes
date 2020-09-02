function adj=comp_adj(tree)

nbnodes=length(tree);
adj=zeros(nbnodes);

for k=1:nbnodes
    for l=k+1:nbnodes
        if(tree(l)==k)
            adj(k,l)=1;
        end
    end
end

adj=triu(adj)+triu(adj,1)';

end