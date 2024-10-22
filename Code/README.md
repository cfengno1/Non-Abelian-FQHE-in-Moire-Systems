Code in Fig.1b calculate the band structure of the noninteracting model for twisted MoTe2 Eq(1):
  The band parameters are in the subroutine make_band(). Output band structure along G-M-K-G is stored in val.dat.

twist_1band.f calculates the energy spectrum, density structure factors, Chern numbers and particle-cut entanglement spectrum of the many-body Hamiltonian Eq.(5):

1.Set parameters "nx0, ny0", to your desired system sizes.  Currently nx0 * ny0=28. Usually one can set the number of sites ns0=nx0*ny0 to be smaller or equal to 32.

2.Set "ne" for different hole band filling numbers.  Currently ne=Ns0/2 for nu=1/2 filling.

3.In the main code,  set parameters "epsilon" for dielectric constant and twist angle "tw_angle". 

4.Set parameters, V1, V2, W1, W2, phi as described in the paper.

5.To calculate energy spectrum and spin structure factor, set nphase=1.

6.To calculate Chern number, set nphase=13 for total number of square meshes 12 x 12 in boundary phase space.

7.The results are in Energy.dat (for energy), Sq.dat (density structure factor), and Chern_num.dat (chern number of the many-body states).

twist_metric.f calculating the quantum metric of the band as shown in Fig.1(c,d).
