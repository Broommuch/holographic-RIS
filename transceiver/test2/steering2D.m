function v = steering2D(theta,phi,Ny,Nz,dy,dz,lambda)
ky = 2*pi/lambda*sin(theta)*cos(phi);
kz = 2*pi/lambda*sin(theta)*sin(phi);
ay = exp(1j*ky*(0:Ny-1).'*dy);
az = exp(1j*kz*(0:Nz-1).'*dz);
v = kron(az,ay);
v = v / norm(v);
end