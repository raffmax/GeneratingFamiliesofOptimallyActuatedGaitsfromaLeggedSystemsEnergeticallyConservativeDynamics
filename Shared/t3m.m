function C_tensor3 = t3m(A_tensor3,B_tensor3)
%t3m 3rd order tensor multiplication
%   A_tensor3 is of dimension a x b x c
%   B_tensor3 is of dimension b x d x c
% https://joss.tcnj.edu/wp-content/uploads/sites/176/2012/04/2009-Miller.pdf

[mA,nA,rA] = size(A_tensor3);
[mB,nB,rB] = size(B_tensor3);

if nA~=mB
   error('The second dimension of A_tensor3 must match the first dimension of B_tensor3!')
end

if rB == 1
   if nB == 1
       % B is a column vector
       C_tensor3 = reshape(reshape(permute(A_tensor3,[1 3 2]),[],nA)*B_tensor3,mA,rA,[]);
   else
       % B is a matrix
       C_tensor3 = permute(reshape(reshape(permute(A_tensor3,[1 3 2]),[],nA)*B_tensor3,mA,rA,[]),[1 3 2]);
   end
else
    warning('needs implementation!')
end

end

