%A is 4th order tensor and B is 2nd order tensor

function R = doubledottf(A,B)

R(1:3,1:3) = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                R(k,l) = R(k,l) + A(i,j)*B(i,j,k,l);
            end
        end
    end
end
end


            
        
