function [b] = make_multi_vector_b(spanA,kerA,betas,onlykernel)
    arguments
        spanA
        kerA
        betas = [0]
        onlykernel = 0
    end
    spanA = spanA';
    kerA = kerA';
    dim_spanA = size(spanA,2);
    dim_kerA = size(kerA,2);
    if norm(kerA) < 1e-15
        dim_kerA = 0;
    end
    b = zeros(dim_spanA+dim_kerA,length(betas));
    randomvectorker = rand(dim_kerA,1);
    randomvectorspan = rand(dim_spanA,1);
    if onlykernel == 1
        randomvectorspan = zeros(dim_spanA,1);
    end
    if norm(kerA) < 1e-10
        for i = 1: size(betas,2)
        b(:,i) = spanA * randomvectorspan;
        end
    else
        norm_ker = norm(kerA*randomvectorker);
        for i = 1: size(betas,2)
            b(:,i) = spanA * randomvectorspan + betas(i) * kerA * randomvectorker/norm_ker;
            %b(:, i) = cos(deltas(i)) * spanA * randomvectorspan/norm(spanA * randomvectorspan) + sin(deltas(i)) * kerA * randomvectorker/norm(kerA * randomvectorker);
        end
    end