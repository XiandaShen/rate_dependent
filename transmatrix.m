function [nP]=transmatrix(psiQ, thetaQ, phiQ,oP)
QQ=[cos(psiQ)*cos(thetaQ)*cos(phiQ)-sin(psiQ)*sin(phiQ)      sin(psiQ)*cos(thetaQ)*cos(phiQ)+cos(psiQ)*sin(phiQ)   -sin(thetaQ)*cos(phiQ)
    -cos(psiQ)*cos(thetaQ)*sin(phiQ)-sin(psiQ)*cos(phiQ)     -sin(psiQ)*cos(thetaQ)*sin(phiQ)+cos(psiQ)*cos(phiQ)  sin(thetaQ)*sin(phiQ)           
    cos(psiQ)*sin(thetaQ)                                  sin(psiQ)*sin(thetaQ)                  cos(thetaQ)];
nP(1:3,1:3,1:3,1:3) = 0;
tQ=QQ';
%tQ=eye(3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for o=1:3
                            for p=1:3
                             nP(i,j,k,l)= nP(i,j,k,l)+  tQ(i,m)*tQ(j,n)*tQ(k,o)*tQ(l,p)*oP(m,n,o,p);
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
                
                
               
        






end
