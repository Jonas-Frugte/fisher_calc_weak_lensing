import data_importer_new as di
import numpy as np
import matplotlib.pyplot as plt

ls = np.arange(2, 2000, 2)

fig, axs = plt.subplots(7, 1, figsize=(6, 21))

pars = [b'H', b'ombh2', b'omch2', b'As', b'mnu', b'ns', b'w0']

for i in range(len(pars)):
    for dd in [-2, -1, 0, 1, 2]:
        axs[i].semilogx(ls, [di.lbs_der_test(l, l, l, b'c', b'c', b'c', 200, pars[i], dd) / di.lbs_der_test(l, l, l, b'c', b'c', b'c', 200, pars[i], 0) for l in ls], label = f'{pars[i]}, {dd}')
        axs[i].legend()

fig.savefig('plots/derivative_accuracy.png')