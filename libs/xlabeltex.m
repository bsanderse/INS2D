function xlabeltex( string, fontsize)
%XLABELTEX xlabel with latex typesetting
% string is the latex string without dollar signs
% fontsize is the size of the latex font

    xlabel(['$' string '$'],'Interpreter','Latex','Fontsize',fontsize);

end

