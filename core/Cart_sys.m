%%%%%% Example RR robot
clc,clear;
syms d1 s2 s3 g dd1 ds2 ds3 ddd1 dds2 dds3 'real' 
syms  a1 a2 a3 I1xx I1yy I1zz I2xx I2yy I2zz I3xx I3yy I3zz 'real'
syms lc1 lc2 lc3 m1 m2 m3 'real'

gc=[0;g;0];
m=[m1;m2;m3];
q=[d1,s2,s3]; dq=[dd1 ds2 ds3]; ddq=[ddd1 dds2 dds3];
% lc1=a1/2;
% lc2=a2/2;
I1=[I1xx 0 0;0 I1yy 0;0 0 I1zz];
I2=[I2xx 0 0;0 I2yy 0;0 0 I2zz];
I3=[I3xx 0 0;0 I3yy 0;0 0 I3zz];

DHpram=[a1,-pi/2,pi/2,d1,0;
        a2,0,s2,0,1;
        a3,0,s3,0,1];

% ICi relative to frame i-1 expressed in frame i 

IC_p0=[lc1,  0, 0;
       lc2,  0, 0;
       lc3,  0, 0];

F1=(m3.*(lc3.*sin(s2 + s3) + a2.*sin(s2)) + lc2.*m2.*sin(s2)).*ds2.^2 + 2.*lc3.*m3.*sin(s2 + s3).*ds2.*ds3 + lc3.*m3.*sin(s2 + s3).*ds3.^2 - dds2.*(m3.*(lc3.*cos(s2 + s3) + a2.*cos(s2)) + lc2.*m2.*cos(s2)) + ddd1.*(m1 + m2 + m3) - dds3.*lc3.*m3.*cos(s2 + s3);
T2=I2zz.*dds2 + I3zz.*dds2 + I3zz.*dds3 + a2.^2.*dds2.*m3 + dds2.*lc2.^2.*m2 + dds2.*lc3.^2.*m3 + dds3.*lc3.^2.*m3 - ddd1.*lc3.*m3.*cos(s2 + s3) + g.*lc3.*m3.*sin(s2 + s3) - a2.*ddd1.*m3.*cos(s2) - ddd1.*lc2.*m2.*cos(s2) + a2.*g.*m3.*sin(s2) + g.*lc2.*m2.*sin(s2) - a2.*ds3.^2.*lc3.*m3.*sin(s3) + 2.*a2.*dds2.*lc3.*m3.*cos(s3) + a2.*dds3.*lc3.*m3.*cos(s3) - 2.*a2.*ds2.*ds3.*lc3.*m3.*sin(s3);
T3=I3zz.*dds2 + I3zz.*dds3 + dds2.*lc3.^2.*m3 + dds3.*lc3.^2.*m3 - ddd1.*lc3.*m3.*cos(s2 + s3) + g.*lc3.*m3.*sin(s2 + s3) + a2.*ds2.^2.*lc3.*m3.*sin(s3) + a2.*dds2.*lc3.*m3.*cos(s3);


I(:,:,1)=I1;
I(:,:,2)=I2;
I(:,:,3)=I3;
[M,V,G]=DH_dyn(DHpram,IC_p0,q,dq,I,m,gc);
tau=simplify(M*[ddd1; dds2;dds3]+V+G);

e1=simplify(tau(1)-F1)
e2=simplify(tau(2)-T2)
e3=simplify(tau(3)-T3)

