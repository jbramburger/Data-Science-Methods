function [U,S,V,threshold,w,sortdog,sortcat] = dc_trainer(dog_wave,cat_wave,feature)
    
    nd = size(dog_wave,2);
    nc = size(cat_wave,2);
    [U,S,V] = svd([dog_wave cat_wave],'econ'); 
    animals = S*V';
    U = U(:,1:feature); % Add this in
    dogs = animals(1:feature,1:nd);
    cats = animals(1:feature,nd+1:nd+nc);
    md = mean(dogs,2);
    mc = mean(cats,2);

    Sw = 0;
    for k=1:nd
        Sw = Sw + (dogs(:,k)-md)*(dogs(:,k)-md)';
    end
    for k=1:nc
        Sw = Sw + (cats(:,k)-mc)*(cats(:,k)-mc)';
    end
    Sb = (md-mc)*(md-mc)';
    
    [V2,D] = eig(Sb,Sw);
    [lambda,ind] = max(abs(diag(D)));
    w = V2(:,ind);
    w = w/norm(w,2);
    vdog = w'*dogs;
    vcat = w'*cats;
    
    if mean(vdog)>mean(vcat)
        w = -w;
        vdog = -vdog;
        vcat = -vcat;
    end
    
    % Don't need plotting here
    sortdog = sort(vdog);
    sortcat = sort(vcat);
    t1 = length(sortdog);
    t2 = 1;
    while sortdog(t1)>sortcat(t2)
    t1 = t1-1;
    t2 = t2+1;
    end
    threshold = (sortdog(t1)+sortcat(t2))/2;

    % We don't need to plot results
end

