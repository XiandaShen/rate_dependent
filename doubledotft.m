%A is 4th order tensor and B is 2nd order tensor

function R = doubledotft(A,B)

R(1:3,1:3) = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                R(i,j) = R(i,j) + A(i,j,k,l)*B(k,l);
            end
        end
    end
end


            
        
