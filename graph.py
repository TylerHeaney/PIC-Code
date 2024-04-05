import numpy as np
import matplotlib.pyplot as plt

x=np.load('pos.npy')
v=np.load('vel.npy')
t=np.load('time.npy')
f=np.load('f.npy')
x=x/100
v=v/100
def calc(p):
    def calc1(p):
        return 1/4/np.pi/8.85418e-12/100 * 1.6022e-19**2 / (p-.1*.00001*400/100)**2 / 1.67e-27
    def calc2(p):
        return -1/4/np.pi/8.85418e-12/100 * 1.6022e-19**2 / (p-.9*.00001*400/100)**2 / 1.67e-27
    return calc1(p)+calc2(p)

fig,ax=plt.subplots()
# ax.plot(t,v)
# ax.plddot(t,x)
a=np.diff(v)/(t[1]-t[0])
# plt.plot(x[:-1],a,color='blue')
# plt.plot(x,[calc(p) for p in x],color='orange')
plt.plot(range(len(f)),f)
plt.show()
# fig,ax=plt.subplots()
# ax.plot(t,x)
# plt.show():w
