function [htx, hty, htx_act, hty_act] = heavy(xcdf, cdf, obj)
    
    basex = 1;
    basey = log(1000);

    ind = find(xcdf>0);

    htx = log(xcdf(ind)) / basex;
    hty = log(1-cdf(ind)) / basey;

    obj.x = xcdf(ind);
    obj = dist_list(obj);

    htx_act = log(obj.x) / basex;
    hty_act = log(1-obj.cdf_y) / basey;
end