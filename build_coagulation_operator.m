function [idxK, idxKp1, wK, wKp1] = build_coagulation_operator(v_node)
% BUILD_COAGULATION_OPERATOR  Fixed-pivot (Hounslow) redistribution weights.
%
% Precomputes the indices/weights used to split the volume of a
% coagulated pair (v_i+v_j) between the two bracketing nodes, per the
% paper's eq.(21). Called once at the start of the simulation since it
% only depends on the (fixed) nodal grid, not on time or temperature.

    n = numel(v_node);
    idxK   = zeros(n*n,1);
    idxKp1 = zeros(n*n,1);
    wK     = zeros(n*n,1);
    wKp1   = zeros(n*n,1);
    c = 0;
    for i = 1:n
        for j = 1:n
            c = c + 1;
            vij = v_node(i) + v_node(j);
            if vij >= v_node(end)
                idxK(c) = n; idxKp1(c) = n; wK(c) = 1; wKp1(c) = 0;
            else
                k = find(v_node <= vij, 1, 'last');
                if k == n
                    idxK(c) = n; idxKp1(c) = n; wK(c) = 1; wKp1(c) = 0;
                else
                    kp1 = k + 1;
                    wk_ = (v_node(kp1) - vij) / (v_node(kp1) - v_node(k));
                    idxK(c)   = k;
                    idxKp1(c) = kp1;
                    wK(c)     = wk_;
                    wKp1(c)   = 1 - wk_;
                end
            end
        end
    end
end
