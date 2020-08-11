%A is 2th order tensor and B is 2nd order tensor

function R = doubledottt(A,B)

R = 0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                R = R + A(i,j)*B(j,i);
            end
        end
    end
end
end


            
        
