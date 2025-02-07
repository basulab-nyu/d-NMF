function [s,c] = graphconncomp(G)
    % [s,c] = graphconncomp(G)
    % Helper function because MATLAB replaced graphconncomp with conncomp
    c = conncomp(graph(G));
    s = max(c);
end