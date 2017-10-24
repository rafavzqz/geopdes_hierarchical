function k_der = conductivity_der_3D(eu)
%CONDUCTIVITY_DER evaluate the temeprature dependent
%conductivity derivative of the material

k_der = cat (1, ...
    reshape (zeros(size(eu)), [1, size(eu)]), ...
    reshape (zeros(size(eu)), [1, size(eu)]), ...
    reshape (zeros(size(eu)), [1, size(eu)]));

for k =1:3
    for i=1:size(eu,1)
        for j=1:size(eu,2)
            if eu(i,j) > 800
                k_der(k,i,j) =  34.0;
            else
                k_der(k,i,j) =  26.7;
            end
        end
    end
end

end

