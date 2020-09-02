function tree=create_tree(nbnodes,type,nbc)

tree=zeros(1,nbnodes);

if isequal(type,'star')
    tree(2:end)=1;
elseif isequal(type,'path')
    tree(2:end)=[1:nbnodes-1];
else
    tree=reg_tree(tree,nbnodes,nbc);
end

end

function tree=reg_tree(tree,nbnodes,nbc)
    for k=1:nbnodes
       tree(k)=floor((k+nbc-2)/nbc); 
    end
end