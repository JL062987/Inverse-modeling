function BD = myBlkdiag(z,n)

    [d1,d2] = size(z);
    zmc = mat2cell(repmat(z,1,n), d1, ones(1,n)*d2);
    BD = blkdiag(zmc{:}); 

end