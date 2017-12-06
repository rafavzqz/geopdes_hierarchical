function [R] = incl(hmsh, hspace)

% Loop over dimensions
for d = 1:hmsh.ndim
    
    n  = length(hspace.space_of_level(1).knots{d}) - hspace.space_of_level(1).degree(d) -1;
    % Make the first level to first level operator
    Rdim{d}{1}{1} = eye(n,n);

    for lev = 1:hspace.nlevels-1
        
        space = hspace.space_of_level(lev);
        
        knots = space.knots{d};
        knots_nextLev = hspace.space_of_level(lev+1).knots{d};
        p = space.degree(d);
        regularity = hspace.regularity(d);

        %Define refined knots and new knotvector
        [~,~,newknots] = kntrefine(knots, 1, p, regularity);

        n = length(knots) - p - 1;
        m = length(newknots);

        % Construct the zero degree global insertion operator
        T{1} = zeros(n+m,n);
        for i = 1:n+m
            for j = 1:n
                if (knots_nextLev(i) >= knots(j) && knots_nextLev(i) < knots(j+1))
                    T{1}(i,j) = 1;
                end
            end
        end

        % Construct the other degree global insertion operators
        for q = 1:p
            T{q+1} = zeros(n+m,n);

            % Add extra zero column to prevent index error
            T{q} = [T{q} zeros(n+m,1)];
            for i = 1:n+m
                for j = 1:n
                    % Split up left and right side for NaN check
                    left  = (knots_nextLev(i+q)-knots(j))/(knots(j+q)-knots(j))*T{q}(i,j);
                    right = (knots(j+q+1)-knots_nextLev(i+q))/(knots(j+q+1)-knots(j+1))*T{q}(i,j+1);
                    if isnan(left) == 1
                        left = 0;
                    elseif isnan(right) == 1
                        right = 0;
                    end
                    T{q+1}(i,j) = left + right;
                end
            end
        end

        % The last global insertion operator is the Refinement Operator of that
        % level
        Rdim{d}{lev}{lev+1}   = T{q+1}';
        % Define L+1th level to L+1th level operator
        Rdim{d}{lev+1}{lev+1} = eye(n+m,n+m); 
    end

    for L2 = 3:hspace.nlevels
        for L1 = 1:L2-2
            Rdim{d}{L1}{L2} = Rdim{d}{L1}{L2-1} * Rdim{d}{L2-1}{L2};
        end
    end
end

R = Rdim{1};
for d = 1:hmsh.ndim-1
    for L2 = 1:hspace.nlevels
        for L1 = 1:L2
            R{L1}{L2} = kron(Rdim{d}{L1}{L2},R{L1}{L2});
        end
    end
end


% Build truncated refinement operator
if hspace.truncated
    R_trunc{1}{1} = R{1}{1};
    for lev = 1:hspace.nlevels
        R_trunc{lev}{lev+1}   = R{lev}{lev+1} * diag(~union(hspace.active{lev+1},hspace.deactivated{lev+1}));
        R_trunc{lev+1}{lev+1} = R{lev+1}{lev+1};
    end   
    for L2 = 3:hspace.nlevels
        for L1 = 1:L2-2
            R_trunc{L1}{L2} = R_trunc{L1}{L2-1} * R_trunc{L1+1}{L2};
        end
    end
    R = R_trunc;
end

end