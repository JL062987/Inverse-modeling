function mysavegcf(Filename,opt)
    
    if nargin == 1
        opt.eps = 0;
    end

    exportgraphics(gcf,strcat(Filename,".png"),'Resolution',600)
    
    if opt.eps == 1
        exportgraphics(gcf,strcat(Filename,".eps"),'BackgroundColor','none','ContentType','vector')
    end

end
