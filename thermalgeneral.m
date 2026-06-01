function [ltype, rtype, ttype, btype] = thermalgeneral(leftside,rightside,topside,bottomside)


switch(leftside)

    case 'Fixed Temperature'

    ltype = 1;
    case 'Convection'

    ltype = 2;
    case 'Insulated'
        
    ltype = 3;


end


switch(rightside)

    case 'Fixed Temperature'
    

    case 'Convection'


    case 'Insulated'

end




switch(topside)

    case 'Fixed Temperature'


    case 'Convection'


    case 'Insulated'

end



switch(bottomside)

    case 'Fixed Temperature'


    case 'Convection'


    case 'Insulated'

end


end