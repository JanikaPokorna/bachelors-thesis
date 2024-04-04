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
    randomvector = rand(dim_kerA,1);
    b = spanA * rand(dim_spanA,1) + delta * kerA * randomvector/norm(kerA * randomvector);
% mozna pak normovat cast ve span A