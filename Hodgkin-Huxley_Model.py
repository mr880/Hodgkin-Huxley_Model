import scipy as sp
import pylab as plt
from scipy.integrate import odeint

# alpha n and Beta n are given in reciprocal msec
# and V is the displacement of the
# membrane potential from its resting value in mv
class HH():
    # ------> Experimental values chosen

    Cm = 1.0
    V_Na = -115
    V_K = 12
    V_l = -10.613  # exact value chosen to make
    # the total ionic current
    # zero at the resting potential
    g_Na = 120
    g_K = 36
    g_l = 0.3

    t = sp.arange(0.0, 450.0, 0.01)
    # ------->


    def alpha_n(self, V):
        return (0.01 * (V + 10) / (sp.exp((V + 10) / 10) - 1))


    def beta_n(self, V):
        return (0.125 * sp.exp(V / 80))


    def alpha_m(self, V):
        return (0.1 * (V + 25) / (sp.exp((V + 25) / 10) - 1))


    def beta_m(self, V):
        return (4 * sp.exp(V / 18))


    def alpha_h(self, V):
        return (0.07 * sp.exp(V / 20))


    def beta_h(self, V):
        return (1 / (sp.exp((V + 30) / 10) + 1))


    def I_K(self, n, V):
        return (self.g_K * n ** 4) * (V - self.V_K)

    def I_Na(self, m, V, h):
        return (self.g_Na * m ** 3) * h * (V - self.V_Na)

    def I_l(self, V):
        return self.g_l * (V - self.V_l)

    def I_inj(self, t):
        return 10 * (t > 100) - 10 * (t > 200) + 35 * (t > 300) - 35 * (t > 400)



    def dvdt(values, t, self):
        V, m, h, n = values

        return (self.I_inj(t) - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_l(V)) / self.Cm

    def main(self):

        values = odeint(self.dvdt, [-65, 0.05, 0.6, 0.32], self.t, args=(self,))
        V = values[:, 0]
        m = values[:, 1]
        h = values[:, 2]
        n = values[:, 3]

        plt.figure()

        plt.title('Hodgkin-Huxley Model')
        plt.plot(self.t, V, 'k')
        plt.ylabel('V (mV)')

        plt.show()

if __name__ == '__main__':
    runner = HH()
    runner.main()