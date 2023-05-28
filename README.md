The MATLAB code to solve the inverse heat conduction problem originates from Ref. Woodbury K A, Najafi H, Monte F D, et al. Inverse Heat Conduction: Ill-Posed Problems[M]. 2nd. Hoboken:John Wiley & Sons, 2023 and is developed to solve the IHCP in three dimension domain.

The methodology is packed into the function fSFSM. Main_Validatio.m illustrates how to use fSFSM function and proves the same ability as the original file for 1-D IHCP.

## 变量说明 Defination of variable

Ns = Nx*Ny, Number of Sensor

Nf = Nx*Ny, Number of Flux

Nt = Number of Time

Qfo is the estimated heat flux

`[Nf, Nt-r+1] = size(Qfo)`

Phi is the heat response to the reference heat flux
`[Ns, Nf, Nt] = size(Phi)`

Y is the measured temperature

`[Ns, Nt+1] = size(Y)`