from matplotlib import pyplot as plt
import numpy as np

# constants
hbar2m = 6.0596
#beta2 = 1.
alpha = 1.
eps = 10.22
sigma = 2.556
beta1 = (8/25/hbar2m * eps*(sigma**12))**0.1

# functions
def u(r, beta1, beta2):
    return (beta1/r)**beta2

def u_p(r, beta1, beta2):
    return -beta2*(beta1**beta2)/r**(beta2+1)

def u_pp(r, beta1, beta2):
    return beta2*(beta2+1)*(beta1**beta2)/r**(beta2+2)

n_div = 100
r = np.linspace(1/n_div, 0.1 + 1/n_div, n_div)
#beta2_arr = np.linspace(-10, 10, 20)
beta2_arr = [5]
for beta2 in beta2_arr:
    H = hbar2m*(u_pp(r, beta1, beta2) + 2*u_p(r, beta1, beta2)/r - r*u_p(r, beta1, beta2)/alpha - u_p(r, beta1, beta2)**2/2) + 4*eps*(sigma**12/r**12 - sigma**6/r**6)
    plt.plot(r, H)
    
plt.legend([f'$\\beta_2 = {beta2}$' for beta2 in beta2_arr], fontsize='small')
#plt.yscale('log')
plt.show()
