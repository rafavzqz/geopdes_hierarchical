function C = capacity(eu, eu_old )
%CAPACITY = volumetric heat capacity for Ti-6Al-4V from K. C. Mills,"Recommended Values of
%Thermophysic Properties for selected commercial alloys", 2002.

C = zeros(size(eu));

for i=1:size(eu,1)
    for j=1:size(eu,2)
        if eu(i,j) <= 950
            density = 4420 - 0.154*(eu(i,j)-25);
            capacity = (0.55 + 0.2*(eu(i,j)/950))*1.0e+03;
            C(i,j) = density * capacity;
        elseif (950 < eu(i,j) && eu(i,j) <= 1650)
            density = 4420 - 0.154*(eu(i,j)-25);
            capacity = (0.64 + 0.09*((eu(i,j)-950)/(1650-950)))*1.0e+03;
            enthalpyJump = 53.0e+03;
            C(i,j) =  density * capacity + ...
            density * enthalpyJump * func_pc_der_mart(eu(i,j), eu_old(i,j));
        else
            density = 3920 - 0.68*(eu(i,j)-1650);
            capacity = (0.83)*1.0e+03;
            enthalpyJump = 286.0e+03;
            C(i,j) =  density * capacity + ...
            density * enthalpyJump * func_pc_der_fus(eu(i,j), eu_old(i,j)); 
        end
    end
end 

end

%Phase Change function and derivative at martensitic-liquid transition
function f_pcDer = func_pc_der_fus(eu, eu_old)

if eu == eu_old
    f_pcDer = 0.0;
else
    f_pcDer = (func_pc_fus(eu) - func_pc_fus(eu_old)) / (eu - eu_old);
end

end

function f_pc = func_pc_fus(eu)

f_pc = zeros(size(eu));

f_pc(eu <= 1650) = 0.0;
f_pc(eu > 1650) = 1.0;

end

%Phase Change function and derivative at solid(alpha)-martensitic(beta) transition
function f_pcDer = func_pc_der_mart(eu, eu_old)

if eu == eu_old
    f_pcDer = 0.0;
else
    f_pcDer = (func_pc_mart(eu) - func_pc_mart(eu_old)) / (eu - eu_old);
end

end

function f_pc = func_pc_mart(eu)

f_pc = zeros(size(eu));

f_pc(eu <= 950) = 0.0;
f_pc(eu > 950) = 1.0;

end

