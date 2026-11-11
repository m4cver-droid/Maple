restart;

with(LinearAlgebra);

L1:=4;

L2:=6;

L3:=8;

X:=[0,4,10,18];

ndof:=8;

ke:=unapply(E*J/l^3*Matrix([[12,6*l,-12,6*l],[6*l,4*l^2,-6*l,2*l^2],[-12,-6*l,12,-6*l],[6*l,2*l^2,-6*l,4*l^2]]),l);

K1:=ke(L1);

K2:=ke(L2);

K3:=ke(L3);

K:=Matrix(ndof,ndof,0):

for e from 1 to 3 do le:=[L1,L2,L3][e]; kl:=ke(le); for i from 1 to 4 do for j from 1 to 4 do K[2*e-2+i,2*e-2+j]:=K[2*e-2+i,2*e-2+j]+kl[i,j]; od; od; od:

K;

F:=Vector([0,-M,0,0,-P-q*L3/2,-q*L3^2/12,-q*L3/2,q*L3^2/12]);

fixed:=[1,7];

free:=[2,3,4,5,6,8];

Kff:=SubMatrix(K,free,free);

Ff:=Vector(6,i->F[free[i]]);

uf:=LinearSolve(Kff,Ff);

U:=Vector(ndof,0):

for i from 1 to 6 do U[free[i]]:=simplify(uf[i]); od:

U;

R:=map(simplify,K . U-F);

RA:=simplify(R[1]);

RB:=simplify(R[7]);

simplify(RA+RB-P-8*q);

simplify(18*RB-10*P-112*q-M);

N1:=1-3*(s/l)^2+2*(s/l)^3;

N2:=s-2*s^2/l+s^3/l^2;

N3:=3*(s/l)^2-2*(s/l)^3;

N4:=-s^2/l+s^3/l^2;

for e from 1 to 3 do le:=[L1,L2,L3][e]; ve[e]:=simplify(subs(l=le,N1)*U[2*e-1]+subs(l=le,N2)*U[2*e]+subs(l=le,N3)*U[2*e+1]+subs(l=le,N4)*U[2*e+2]); theta[e]:=simplify(diff(ve[e],s)); moment[e]:=simplify(E*J*diff(ve[e],s$2)); shear[e]:=simplify(diff(moment[e],s)); od:

seq(ve[e],e=1..3);

seq(theta[e],e=1..3);

seq(moment[e],e=1..3);

seq(shear[e],e=1..3);
