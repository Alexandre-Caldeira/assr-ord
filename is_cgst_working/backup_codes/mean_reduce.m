function red = mean_reduce(vec,dim)
    red = squeeze(mean(vec,dim,'omitnan'));
end