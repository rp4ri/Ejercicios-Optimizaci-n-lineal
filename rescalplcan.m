function [ch,Ah,bh,uh,r0,r,s,s0] = rescalplcan(c,A,b,u)
[m,n] =  size(A); ind1n=1:n; ind1n=ind1n(:); ind1m=1:m; ind1m=ind1m(:);
c=c(:); b=b(:); u=u(:); can = u < 1.0e+20;
uh = 1.0e+20*ones(n,1);
r0 = max(abs(c)); 
Rx = max(abs(A),[],2); Rx=Rx(:); 
r = max([Rx abs(b)],[],2); r(r==0)=1; r=1./r;
R = spconvert([ind1m ind1m r]);
rA = sparse(R*A); rb=r.*b; r0c=c/r0;
sx = max(abs(rA),[],1); s = max([abs(r0c)'; sx],[],1); s(s==0)=1;
sc = s(can); sc=sc(:);
s = 1./s; S=spconvert([ind1n ind1n s']);
s0 = max(abs(rb));
if s0 == 0, s0 = 1; end
ch = r0c'.*s; ch = ch(:);
Ah = sparse(rA*S); 
bh = rb/s0;
uh(can) = sc.*u(can)/s0;

% s0 =max(abs(b));
% if s0 == 0, s0 = 1; end
% bs = b/s0; bs = bs(:);
% Sx = max(abs(A),[],1); S = max([abs(c)'; Sx],[],1);
% cS = c./S';
% r0 = max(abs(cS));
% ch = cS/r0; ch=ch(:);
% As = sparse(A./S);
% Rx=max(abs(As),[],2); Rx=Rx(:); 
% R = max([Rx abs(bs)],[],2); R=R(:);
% Ah = sparse(As./R);
% bh = bs./R; bh = bh(:); sc = S(can); sc=sc(:);
% uh(can) = sc.*u(can)/s0; uh = uh(:);