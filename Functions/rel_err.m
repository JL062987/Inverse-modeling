function z = rel_err(x,y)

    z = norm(vec(x-y))/norm(vec(y));

end