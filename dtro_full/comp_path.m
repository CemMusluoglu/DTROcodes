function path_tree=comp_path(tree,leafs)

nb_leafs=length(leafs);
path_tree=[];
paths={};

for k=1:nb_leafs
    current=leafs(k);
    path_curr=[current];
    while current~=1
        path_curr=[path_curr,tree(current)];
        current=tree(current);
    end
    paths{k}=flip(path_curr);
end


for k=1:nb_leafs
    
    path_c=paths{k};
    if(k==1)
        path_tree=[path_tree,path_c];
    end
    
    if(k<nb_leafs)
        path_n=paths{k+1};
        common_p=intersect(path_c,path_n);
        transit_node=common_p(end);
        transit_path=path_c(end-1:-1:find(path_c==transit_node,1));
        transit_path=[transit_path,path_n(find(path_n==transit_node)+1:end)];
        path_tree=[path_tree,transit_path];
    else
        transit_path=path_c(end-1:-1:1);
        path_tree=[path_tree,transit_path];
    end
    
    
end

path_tree=path_tree(1:end-1);
path_tree=unique(path_tree,'stable');


end