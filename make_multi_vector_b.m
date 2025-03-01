function [b] = make_multi_vector_b(spanA,kerA,deltas)
    arguments
        spanA
        kerA
        deltas = [0]
    end
    spanA = spanA';
    kerA = kerA';
    dim_spanA = size(spanA,2);
    dim_kerA = size(kerA,2);
    b = zeros(dim_spanA+dim_kerA,length(deltas));
    randomvectorker = rand(dim_kerA,1);
    randomvectorspan = rand(dim_spanA,1);
    for i = 1: size(deltas,2)
        b(:,i) = spanA * randomvectorspan + deltas(i) * kerA * randomvectorker/norm(kerA * randomvectorker);
        %b(:, i) = cos(deltas(i)) * spanA * randomvectorspan/norm(spanA * randomvectorspan) + sin(deltas(i)) * kerA * randomvectorker/norm(kerA * randomvectorker);
    end