function [c,A,b,lo,hi,row_zeros,col_zeros] = presolving(c,A,b,lo,hi)
eps_tol=1.0e-16; b=b(:); c=c(:);
rA = sum(abs(A) > eps_tol,2); 
rb = abs(b) > eps_tol;
rowz=rA==0 & rb==0; row_zeros=full(sum(rowz));
if row_zeros > 0
    A = A(~rowz,:);
    b = b(~rowz);
end
cA = sum(abs(A) > eps_tol,1); cA = cA(:);
cc = abs(c) > eps_tol;
colz=cA==0 & cc==0; col_zeros=full(sum(colz));
if col_zeros > 0
    A = A(:,~colz);
    c = c(~colz);
    if ~isempty(lo), lo = lo(~colz); end
    if ~isempty(hi), hi = hi(~colz); end
end
