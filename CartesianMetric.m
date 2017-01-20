function h=CartesianMetric(X)
%
% h=CartesianMetric(X)
%
% Computes cartesian components of the radiation field 
% (far field metric perturbation h^{jk}) given a worldline 
% X.  The worldline structure X must be produced by KludgedXofl.
%
% Output structure is:
%
%    h.quad.jk = [normalized quadrupole term]   
%    h.oct.ijk = [normalized un-projected part of octupole term]
%    
% for i,j,k in {x, y, z} so tha the total quadrupole-octupole metric 
% perturbation is
%
%  h.jk = (mu/r)[h.quad.jk + n_i h.quad.ijk]
%
% where n_i is the unit 3-vector pointing from the source to the observer.
% This routine doesn't know where the observer is, so it leaves the
% octupole term in an un-projected form.
%
% Use this as input to HPHX in order to get the waveform for a 
% specific observer.
%
% See also KLUDGEDINSPIRAL KLUDGEDXOFL HPHX
%
% Steve Drasco
%

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% quadrupole moment I^{jk} 
Ixx = X.x.^2;
Iyy = X.y.^2;
Izz = X.z.^2;
Ixy = X.x .* X.y;
Ixz = X.x .* X.z;
Iyz = X.y .* X.z;

% first t-derivatives of I^{jk} 
dIxx = 2 * X.x .* X.dx;
dIyy = 2 * X.y .* X.dy;
dIzz = 2 * X.z .* X.dz;
dIxy = (X.dx .* X.y) + (X.x .* X.dy);
dIxz = (X.dx .* X.z) + (X.x .* X.dz);
dIyz = (X.dy .* X.z) + (X.y .* X.dz);

% second t-derivatives of I^{jk} 
ddIxx = 2*(X.dx.^2 + X.x .* X.ddx);
ddIyy = 2*(X.dy.^2 + X.y .* X.ddy);
ddIzz = 2*(X.dz.^2 + X.z .* X.ddz);
ddIxy = (X.ddx .* X.y) + (2 * X.dx .* X.dy) + (X.x .* X.ddy);
ddIxz = (X.ddx .* X.z) + (2 * X.dx .* X.dz) + (X.x .* X.ddz);
ddIyz = (X.ddy .* X.z) + (2 * X.dy .* X.dz) + (X.y .* X.ddz);

% third t-derivatives of I^{jk} 
dddIxx = (2 * X.dddx .* X.x) + (6 * X.ddx .* X.dx);
dddIyy = (2 * X.dddy .* X.y) + (6 * X.ddy .* X.dy);
dddIzz = (2 * X.dddz .* X.z) + (6 * X.ddz .* X.dz);
dddIxy = (X.dddx .* X.y) + (X.x .* X.dddy) ...
    + (3 * X.ddx .* X.dy) + (3 * X.dx .* X.ddy);
dddIxz = (X.dddx .* X.z) + (X.x .* X.dddz) ...
    + (3 * X.ddx .* X.dz) + (3 * X.dx .* X.ddz);
dddIyz = (X.dddy .* X.z) + (X.y .* X.dddz) ...
    + (3 * X.ddy .* X.dz) + (3 * X.dy .* X.ddz);

% third t-derivative of mass octupole moment M^{ijk}
dddMxxx = (X.dddx .* Ixx) + (X.x .* dddIxx) ...
    + (3 * X.ddx .* dIxx) + (3 * X.dx .* ddIxx);
dddMyxx = (X.dddy .* Ixx) + (X.y .* dddIxx) ...
    + (3 * X.ddy .* dIxx) + (3 * X.dy .* ddIxx);
dddMzxx = (X.dddz .* Ixx) + (X.z .* dddIxx) ...
    + (3 * X.ddz .* dIxx) + (3 * X.dz .* ddIxx);
clear dddIxx;
dddMxyy = (X.dddx .* Iyy) + (X.x .* dddIyy) ...
    + (3 * X.ddx .* dIyy) + (3 * X.dx .* ddIyy);
dddMyyy = (X.dddy .* Iyy) + (X.y .* dddIyy) ...
    + (3 * X.ddy .* dIyy) + (3 * X.dy .* ddIyy);
dddMzyy = (X.dddz .* Iyy) + (X.z .* dddIyy) ...
    + (3 * X.ddz .* dIyy) + (3 * X.dz .* ddIyy);
clear dddIyy;
dddMxzz = (X.dddx .* Izz) + (X.x .* dddIzz) ...
    + (3 * X.ddx .* dIzz) + (3 * X.dx .* ddIzz);
dddMyzz = (X.dddy .* Izz) + (X.y .* dddIzz) ...
    + (3 * X.ddy .* dIzz) + (3 * X.dy .* ddIzz);
dddMzzz = (X.dddz .* Izz) + (X.z .* dddIzz) ...
    + (3 * X.ddz .* dIzz) + (3 * X.dz .* ddIzz);
clear dddIzz;
dddMxxy = (X.dddx .* Ixy) + (X.x .* dddIxy) ...
    + (3 * X.ddx .* dIxy) + (3 * X.dx .* ddIxy);
dddMyxy = (X.dddy .* Ixy) + (X.y .* dddIxy) ...
    + (3 * X.ddy .* dIxy) + (3 * X.dy .* ddIxy);
dddMzxy = (X.dddz .* Ixy) + (X.z .* dddIxy) ...
    + (3 * X.ddz .* dIxy) + (3 * X.dz .* ddIxy);
clear dddIxy;
dddMxxz = (X.dddx .* Ixz) + (X.x .* dddIxz) ...
    + (3 * X.ddx .* dIxz) + (3 * X.dx .* ddIxz);
dddMyxz = (X.dddy .* Ixz) + (X.y .* dddIxz) ...
    + (3 * X.ddy .* dIxz) + (3 * X.dy .* ddIxz);
dddMzxz = (X.dddz .* Ixz) + (X.z .* dddIxz) ...
    + (3 * X.ddz .* dIxz) + (3 * X.dz .* ddIxz);
clear dddIxz;
dddMxyz = (X.dddx .* Iyz) + (X.x .* dddIyz) ...
    + (3 * X.ddx .* dIyz) + (3 * X.dx .* ddIyz);
dddMyyz = (X.dddy .* Iyz) + (X.y .* dddIyz) ...
    + (3 * X.ddy .* dIyz) + (3 * X.dy .* ddIyz);
dddMzyz = (X.dddz .* Iyz) + (X.z .* dddIyz) ...
    + (3 * X.ddz .* dIyz) + (3 * X.dz .* ddIyz);
clear dddIyz;

% second t-derivatives of current quadrupole moment S^{ijk}
ddSxxx = (X.dddx .* Ixx) + (2 * X.ddx .* dIxx) + (X.dx .* ddIxx);
ddSyxx = (X.dddy .* Ixx) + (2 * X.ddy .* dIxx) + (X.dy .* ddIxx);
ddSzxx = (X.dddz .* Ixx) + (2 * X.ddz .* dIxx) + (X.dz .* ddIxx);
clear Ixx dIxx;
ddSxyy = (X.dddx .* Iyy) + (2 * X.ddx .* dIyy) + (X.dx .* ddIyy);
ddSyyy = (X.dddy .* Iyy) + (2 * X.ddy .* dIyy) + (X.dy .* ddIyy);
ddSzyy = (X.dddz .* Iyy) + (2 * X.ddz .* dIyy) + (X.dz .* ddIyy);
clear Iyy dIyy;
ddSxzz = (X.dddx .* Izz) + (2 * X.ddx .* dIzz) + (X.dx .* ddIzz);
ddSyzz = (X.dddy .* Izz) + (2 * X.ddy .* dIzz) + (X.dy .* ddIzz);
ddSzzz = (X.dddz .* Izz) + (2 * X.ddz .* dIzz) + (X.dz .* ddIzz);
clear Izz dIyy;
ddSxxy = (X.dddx .* Ixy) + (2 * X.ddx .* dIxy) + (X.dx .* ddIxy);
ddSyxy = (X.dddy .* Ixy) + (2 * X.ddy .* dIxy) + (X.dy .* ddIxy);
ddSzxy = (X.dddz .* Ixy) + (2 * X.ddz .* dIxy) + (X.dz .* ddIxy);
clear Ixy dIxy;
ddSxxz = (X.dddx .* Ixz) + (2 * X.ddx .* dIxz) + (X.dx .* ddIxz);
ddSyxz = (X.dddy .* Ixz) + (2 * X.ddy .* dIxz) + (X.dy .* ddIxz);
ddSzxz = (X.dddz .* Ixz) + (2 * X.ddz .* dIxz) + (X.dz .* ddIxz);
clear Ixz dIxz;
ddSxyz = (X.dddx .* Iyz) + (2 * X.ddx .* dIyz) + (X.dx .* ddIyz);
ddSyyz = (X.dddy .* Iyz) + (2 * X.ddy .* dIyz) + (X.dy .* ddIyz);
ddSzyz = (X.dddz .* Iyz) + (2 * X.ddz .* dIyz) + (X.dz .* ddIyz);
clear Iyz dIyz;

% octupole term, clearing things as we go
h.oct.xxx = -2*ddSxxx + dddMxxx;
clear ddSxxx dddMxxx;
h.oct.xyy = -2*ddSxyy + dddMxyy;
clear ddSxyy dddMxyy;
h.oct.xzz = -2*ddSxzz + dddMxzz;
clear ddSxzz dddMxzz;
h.oct.xxy = -2*ddSxxy + dddMxxy;
clear ddSxxy dddMxxy;
h.oct.xxz = -2*ddSxxz + dddMxxz;
clear ddSxxz dddMxxz;
h.oct.xyz = -2*ddSxyz + dddMxyz;
clear ddSxyz dddMxyz;
h.oct.yxx = -2*ddSyxx + dddMyxx;
clear ddSyxx dddMyxx;
h.oct.yyy = -2*ddSyyy + dddMyyy;
clear ddSyyy dddMyyy;
h.oct.yzz = -2*ddSyzz + dddMyzz;
clear ddSyzz dddMyzz;
h.oct.yxy = -2*ddSyxy + dddMyxy;
clear ddSyxy dddMyxy;
h.oct.yxz = -2*ddSyxz + dddMyxz;
clear ddSyxz dddMyxz;
h.oct.yyz = -2*ddSyyz + dddMyyz;
clear ddSyyz dddMyyz;
h.oct.zxx = -2*ddSzxx + dddMzxx;
clear ddSzxx dddMzxx;
h.oct.zyy = -2*ddSzyy + dddMzyy;
clear ddSzyy dddMzyy;
h.oct.zzz = -2*ddSzzz + dddMzzz;
clear ddSzzz dddMzzz;
h.oct.zxy = -2*ddSzxy + dddMzxy;
clear ddSzxy dddMzxy;
h.oct.zxz = -2*ddSzxz + dddMzxz;
clear ddSzxz dddMzxz;
h.oct.zyz = -2*ddSzyz + dddMzyz;
clear ddSzyz dddMzyz;

% quadrupole term
h.quad.xx = 2 * ddIxx;
clear ddIxx;
h.quad.yy = 2 * ddIyy;
clear ddIyy;
h.quad.zz = 2 * ddIzz;
clear ddIzz;
h.quad.xy = 2 * ddIxy;
clear ddIxy;
h.quad.xz = 2 * ddIxz;
clear ddIxz;
h.quad.yz = 2 * ddIyz;
clear ddIyz;

% log the computational cost of this job
h.CPUsec = cputime - InitialCPUTime;
