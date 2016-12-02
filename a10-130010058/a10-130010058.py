from pysph.solver.utils import load
import numpy as np
from matplotlib import pyplot as plt

# There are total twenty seven cases, so forming a data matrix of 27 elements.


def P_exact(x, y):
    decay_rate = -8.0*np.pi*np.pi/100.0
    t = 2.0
    U = 1.0
    L = 1.0
    rho0 = 1.0
    c0 = 10.0*U
    p0 = c0*c0*rho0
    factor = U*np.exp(decay_rate*t)
    u = factor*(-np.cos(2.0*np.pi*x)*np.sin(2.0*np.pi*y))
    v = factor*(np.sin(2.0*np.pi*x)*np.cos(2.0*np.pi*y))
    p_e = factor*factor*(-0.25*(np.cos(4.0*np.pi*x) + np.cos(4.0*np.pi*y)))
    return u, v, p0+p_e

paths = ['./Cubicspline/case1/taylor_green_4400.npz',
         './Cubicspline/case2/taylor_green_8800.npz',
         './Cubicspline/case3/taylor_green_17600.npz',
         './Cubicspline/case4/taylor_green_2200.npz',
         './Cubicspline/case5/taylor_green_4400.npz',
         './Cubicspline/case6/taylor_green_8800.npz',
         './Cubicspline/case7/taylor_green_1100.npz',
         './Cubicspline/case8/taylor_green_2200.npz',
         './Cubicspline/case9/taylor_green_4400.npz',
         './gaussian/case1/taylor_green_4400.npz',
         './gaussian/case2/taylor_green_8800.npz',
         './gaussian/case3/taylor_green_17600.npz',
         './gaussian/case4/taylor_green_2200.npz',
         './gaussian/case5/taylor_green_4400.npz',
         './gaussian/case6/taylor_green_8800.npz',
         './gaussian/case7/taylor_green_1100.npz',
         './gaussian/case8/taylor_green_2200.npz',
         './gaussian/case9/taylor_green_4400.npz',
         './quintic/case1/taylor_green_4400.npz',
         './quintic/case2/taylor_green_8800.npz',
         './quintic/case3/taylor_green_17600.npz',
         './quintic/case4/taylor_green_2200.npz',
         './quintic/case5/taylor_green_4400.npz',
         './quintic/case6/taylor_green_8800.npz',
         './quintic/case7/taylor_green_1100.npz',
         './quintic/case8/taylor_green_2200.npz',
         './quintic/case9/taylor_green_4400.npz',
         './WendlandQuintic/case1/taylor_green_4400.npz',
         './WendlandQuintic/case2/taylor_green_8800.npz',
         './WendlandQuintic/case3/taylor_green_17600.npz',
         './WendlandQuintic/case4/taylor_green_2200.npz',
         './WendlandQuintic/case5/taylor_green_4400.npz',
         './WendlandQuintic/case6/taylor_green_8800.npz',
         './WendlandQuintic/case7/taylor_green_1100.npz',
         './WendlandQuintic/case8/taylor_green_2200.npz',
         './WendlandQuintic/case9/taylor_green_4400.npz']

Err_p = np.zeros(len(paths))
Err_vl1 = np.zeros(len(paths))
Err_vlinf = np.zeros(len(paths))

for i, path in enumerate(paths):
    data = load(path)
    particle_arrays = data['arrays']
    fluid = particle_arrays['fluid']
    p = fluid.p
    x = fluid.x
    y = fluid.y
    u = fluid.u
    v = fluid.v
    v_mag = np.sqrt(u*u + v*v)
    u_e, v_e, p_e = P_exact(x, y)
    v_mag_e = np.sqrt(u_e*u_e + v_e*v_e)
    Err_vlinf[i] = abs(v_mag_e.max() - v_mag.max())/v_mag_e.max()
    Err_vl1[i] = np.average(np.abs(v_mag-v_mag_e))/np.average(np.abs(v_mag_e))
    Err_p[i] = np.average(abs(p_e - p))/p_e.max()

if np.argmin(Err_p) < 9:
    stringp = 'Cubic Spline'
elif np.argmin(Err_p) >= 9 and np.argmin(Err_p) < 18:
    stringp = 'Gaussian Spline'
elif np.argmin(Err_p) >= 18 and np.argmin(Err_p) < 27:
    stringp = 'Quintic Spline'
else:
    stringp = 'WendlandQuintic'
print "The case with least L_1 error in  pressure is case " + str(np.argmin(Err_p) % 9 + 1) + ' of\
 ' + stringp
if np.argmin(Err_vl1) < 9:
    stringvl1 = 'Cubic Spline'
elif np.argmin(Err_vl1) >= 9 and np.argmin(Err_vl1) < 18:
    stringvl1 = 'Gaussian Spline'
elif np.argmin(Err_vl1) >= 18 and np.argmin(Err_vl1) < 27:
    stringvl1 = 'Quintic Spline'
else:
    stringvl1 = 'WendlandQuintic'
print "The case with least L_1 error in velocity magnitude is case " + str(np.argmin(Err_vl1) % 9 + 1) + ' of\
 ' + stringvl1
if np.argmin(Err_vlinf) < 9:
    stringvlinf = 'Cubic Spline'
elif np.argmin(Err_vlinf) >= 9 and np.argmin(Err_vlinf) < 18:
    stringvlinf = 'Gaussian Spline'
elif np.argmin(Err_vlinf) >= 18 and np.argmin(Err_vlinf) < 27:
    stringvlinf = 'Quintic Spline'
else:
    stringvlinf = 'WendlandQuintic'
print "The case with least L_inf error in velocity magnitude is case " + str(np.argmin(Err_vlinf) % 9 + 1) + ' of\
 ' + stringvlinf

data1 = load('./Cubicspline/case2/taylor_green_8800.npz')
data2 = load('./gaussian/case2/taylor_green_8800.npz')
data3 = load('./quintic/case2/taylor_green_8800.npz')
data4 = load('./WendlandQuintic/case2/taylor_green_8800.npz')
Err = np.zeros(4)
for i, data in enumerate([data1, data2, data3, data4]):
    particle_arrays = data['arrays']
    fluid = particle_arrays['fluid']
    p = fluid.p
    x = fluid.x
    y = fluid.y
    u_e, v_e, p_e = P_exact(x, y)
    Err[i] = np.average(abs(p_e - p))/p_e.max()
plt.bar([1, 2, 3, 4], [Err[0], Err[1], Err[2], Err[3]], align='center')
plt.xticks([1, 2, 3, 4], ["Cubic", "Gaussian", "Quintic", "WendtlQuintic"])
plt.title('L_1 error in Pressure for hdx = 0.5, nx = 50')
plt.savefig('perror.png')
plt.close()
data1 = load('./Cubicspline/case9/taylor_green_4400.npz')
data2 = load('./gaussian/case9/taylor_green_4400.npz')
data3 = load('./quintic/case9/taylor_green_4400.npz')
data4 = load('./WendlandQuintic/case9/taylor_green_4400.npz')
Errv1 = np.zeros(4)
Errv2 = np.zeros(4)
for i, data in enumerate([data1, data2, data3, data4]):
    particle_arrays = data['arrays']
    fluid = particle_arrays['fluid']
    x = fluid.x
    y = fluid.y
    u = fluid.u
    v = fluid.v
    v_mag = np.sqrt(u*u + v*v)
    u_e, v_e, p_e = P_exact(x, y)
    v_mag_e = np.sqrt(u_e*u_e + v_e*v_e)
    Errv1[i] = abs(v_mag_e.max() - v_mag.max())/v_mag_e.max()
    Errv2[i] = np.average(np.abs(v_mag-v_mag_e))/np.average(np.abs(v_mag_e))
plt.bar([1, 2, 3, 4], [Errv1[0], Errv1[1], Errv1[2], Errv1[3]], align='center')
plt.xticks([1, 2, 3, 4], ["Cubic", "Gaussian", "Quintic", "WendtlQuintic"])
plt.title('L_inf error in Velocity magnitude for hdx = 2.0, nx = 100')
plt.savefig('vinferror.png')
plt.close()
plt.bar([1, 2, 3, 4], [Errv2[0], Errv2[1], Errv2[2], Errv2[3]], align='center')
plt.xticks([1, 2, 3, 4], ["Cubic", "Gaussian", "Quintic", "WendtlQuintic"])
plt.title('L_1 error in Velocity magnitude for hdx = 2.0, nx = 100')
plt.savefig('v1error.png')
plt.close()
