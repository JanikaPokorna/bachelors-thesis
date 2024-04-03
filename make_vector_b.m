function [b] = make_vector_b(spanA,kerA,delta)
    arguments
        spanA
        kerA
        delta = 1e-3
    end
    spanA = spanA';
    kerA = kerA';
    dim_spanA = size(spanA,2);
    dim_kerA = size(kerA,2);
    b = spanA * rand(dim_spanA,1) + delta * kerA * rand(dim_kerA,1)/norm(kerA * rand(dim_kerA,1));
