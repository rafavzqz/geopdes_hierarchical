function k = conductivity(eu)
%CONDUCTIVITY evaluate the temeprature dependent 
%conductivity of the material

k = zeros(size(eu));

for i=1:size(eu,1)
    for j=1:size(eu,2)
        if eu(i,j) <= 950
            k(i,j) =  7.5 + 15.0 * (eu(i,j)/950);
        elseif (950 < eu(i,j) && eu(i,j) <= 1650)
            k(i,j) =  12.5 + (27.5-12.5) * ((eu(i,j)-950)/(1650-950));
        else
            k(i,j) =  34.0;
        end
    end
end

% for i=1:size(eu,1)
%     for j=1:size(eu,2)
%         if eu(i,j) >= 800
%             k(i,j) =  34.0 * (eu(i,j)/800);
%         else
%             k(i,j) =  34; % 26.7 * ((800-eu(i,j))/800);
%         end
%     end
% end


end

