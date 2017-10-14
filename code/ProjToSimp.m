function proj_vec = ProjToSimp(vec)
    len_vec = length(vec);
    sort_vec = transpose(sort(vec, 'descend'));
    K = fliplr(find(sort_vec - (cumsum(sort_vec) - 1)./(1 : len_vec) > 0));
    K = K(1);
    tau = (sum(sort_vec(1 : K)) - 1)/K;
    proj_vec = vec - tau;
    index_neg = find(proj_vec < 0);
    proj_vec(index_neg) = 0;
end