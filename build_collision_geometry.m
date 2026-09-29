function G = build_collision_geometry(v_node)
% BUILD_COLLISION_GEOMETRY  Volume-dependent part of the Brownian
% collision frequency, eq.(16): beta(i,j) = const * sqrt(T) * G(i,j).
%
% Precomputed once since it only depends on the (fixed) nodal grid; the
% temperature dependence is applied separately at each time step.

    n = numel(v_node);
    G = zeros(n,n);
    for i = 1:n
        for j = 1:n
            G(i,j) = sqrt(1/v_node(i) + 1/v_node(j)) * ...
                     (v_node(i)^(1/3) + v_node(j)^(1/3))^2;
        end
    end
end
