function z = rel_err(x,y)

z = norm(vec(x)-vec(y))/norm(vec(y));

end