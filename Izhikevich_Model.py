from numpy import *
from pylab import *

# Eugene Izhikevich model

# returns dv/dt = 0.04v^2 + 5v + 140 - u + I
def dvdt(v, u, I):
    dvdt = 0.04 * v**2 + 5 * v + 140 - u + I
    return dvdt

#setting our four parameters that are used to alter cell's behavior
a = 0.02
b = 0.2
c = -65
d = 6

# returns du/dt = a(bv - u)
def dudt(v, u):
    dudt = a * (b * v - u)
    return dudt

# model time length (msec)
time_length = 100

# time step
dt = 0.25

# time
t = arange(0, time_length, dt)

# create two arrays of zeros of size length time
V = zeros(len(t))
U = zeros(len(t))

#initialize V and U
V[0] = -70
U[0] = 0.2 * V[0]


for i in range(1, 380):
    if i < 350:
        I = 14
    else:
        I = 0

    V[i] = V[i - 1] + dt * dvdt(V[i - 1], U[i - 1], I)
    U[i] = U[i - 1] + dt * dudt(V[i - 1], U[i - 1])

    # When v reaches 30 mV, the cell fires,
    # and then v is reset to c and u is increased by d
    if (V[i] >= 30):
        V[i] = c
        U[i] = U[i] + d
    i = i + 1

plot(t, V)
title('Izhikevich Model')
ylabel('Membrane Potential (V)')
xlabel('Time (msec)')
show()