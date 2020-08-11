%get the inverse of 4th order tensor F
%4th order tensor has 18 nonzero components, symmetry & transverse isotropy

function Finverse = inversegeneral(F)

Fmatrix = [F(1,1,1,1) F(1,1,2,2) F(1,1,3,3);...
           F(2,2,1,1) F(2,2,2,2) F(2,2,3,3);...
           F(3,3,1,1) F(3,3,2,2) F(3,3,3,3)];
       
W = inv(Fmatrix);

Finverse(1:3,1:3,1:3,1:3) = 0;

Finverse(1,1,1,1) = W(1,1);
Finverse(1,1,2,2) = W(1,2);
Finverse(1,1,3,3) = W(1,3);
Finverse(2,2,1,1) = W(2,1);
Finverse(2,2,2,2) = W(2,2);
Finverse(2,2,3,3) = W(2,3);
Finverse(3,3,1,1) = W(3,1);
Finverse(3,3,2,2) = W(3,2);
Finverse(3,3,3,3) = W(3,3);

Finverse(1,2,1,2) = 1/4/F(1,2,1,2);
Finverse(2,3,2,3) = 1/4/F(2,3,2,3);
Finverse(1,3,1,3) = 1/4/F(1,3,1,3);

Finverse(1,2,2,1) = Finverse(1,2,1,2);
Finverse(2,1,1,2) = Finverse(1,2,1,2);
Finverse(2,1,2,1) = Finverse(1,2,1,2);

Finverse(1,3,3,1) = Finverse(1,3,1,3);
Finverse(3,1,1,3) = Finverse(1,3,1,3);
Finverse(3,1,3,1) = Finverse(1,3,1,3);

Finverse(3,2,3,2) = Finverse(2,3,2,3);
Finverse(3,2,2,3) = Finverse(2,3,2,3);
Finverse(2,3,3,2) = Finverse(2,3,2,3);





