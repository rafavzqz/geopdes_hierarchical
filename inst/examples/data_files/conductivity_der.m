function k_der = conductivity_der(eu)
%CONDUCTIVITY_DER evaluate the temprature dependent
%conductivity derivative of the material

k_der = cat (1, ...
    reshape (zeros(size(eu)), [1, size(eu)]), ...
    reshape (zeros(size(eu)), [1, size(eu)]));

for k =1:2
    for i=1:size(eu,1)
        for j=1:size(eu,2)
            if eu(i,j) > 800
                k_der(k,i,j) =  34.0/800;
            else
                k_der(k,i,j) =  0.0; %-26.7/800;
            end
        end
    end
end

end

