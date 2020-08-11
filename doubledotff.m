%both A and B are 4th order tensors

function R = doubledotff(A,B)

R(1:3,1:3,1:3,1:3) = 0;

for i = 1:3
    for j = 1:3
        for m = 1:3
            for n = 1:3
%                 R(i,j,m,n)=0;
                for k = 1:3
                    for l = 1:3
                        R(i,j,m,n) = R(i,j,m,n) + A(i,j,k,l)*B(k,l,m,n);
                    end
                end
            end
        end
    end
end

