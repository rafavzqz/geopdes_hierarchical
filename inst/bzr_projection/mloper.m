function [ml] = mloper(activeFuns,lev,R)

ml = cell(1, lev);
for l = 1:lev

    ml{l} = [];
    
    for k = 1:l
               
        aux = R{k}{l}(activeFuns{k});
        ml{l} = [ml{l}; aux];
    end

end

end