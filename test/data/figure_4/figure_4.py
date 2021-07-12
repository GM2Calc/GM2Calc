import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import scipy.interpolate
from matplotlib.backends.backend_pdf import PdfPages

def get_grid(filename):
    x, y, z = np.genfromtxt(filename, unpack=True, dtype=float)
    x = x[~np.isnan(z)]
    y = y[~np.isnan(z)]
    z = z[~np.isnan(z)]
    z = z*1e10

    N = 200
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

    return xi, yi, zi

with PdfPages(r'figure_4.pdf') as pdf:
    xi, yi, zi = get_grid(r'figure_4a.txt')
    fig = plt.figure(figsize=(5,5))
    cs = plt.contour(xi, yi, zi, levels = [-2, -1.5, -1, -0.5, 0])
    plt.clabel(cs, cs.levels, inline=True, fmt="%2.2g", fontsize=10)
    plt.title(r"$M_A = 10$ GeV")
    plt.xlabel(r"$M_{H^+}$ / GeV")
    plt.ylabel(r"$M_H$ / GeV")
    pdf.savefig()
    plt.close()

    xi, yi, zi = get_grid(r'figure_4b.txt')
    fig = plt.figure(figsize=(5,5))
    cs = plt.contour(xi, yi, zi, levels = [-2, -1.5, -1, -0.5, 0])
    plt.clabel(cs, cs.levels, inline=True, fmt="%2.2g", fontsize=10)
    plt.title(r"$M_A = 50$ GeV")
    plt.xlabel(r"$M_{H^+}$ / GeV")
    plt.ylabel(r"$M_H$ / GeV")
    pdf.savefig()
    plt.close()

    xi, yi, zi = get_grid(r'figure_4c.txt')
    fig = plt.figure(figsize=(5,5))
    cs = plt.contour(xi, yi, zi, levels = [-1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25])
    plt.clabel(cs, cs.levels, inline=True, fmt="%2.2g", fontsize=10)
    plt.title(r"$M_A = 100$ GeV")
    plt.xlabel(r"$M_{H^+}$ / GeV")
    plt.ylabel(r"$M_H$ / GeV")
    pdf.savefig()
    plt.close()

    xi, yi, zi = get_grid(r'figure_4d.txt')
    fig = plt.figure(figsize=(5,5))
    cs = plt.contour(xi, yi, zi, levels = [-1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2])
    plt.clabel(cs, cs.levels, inline=True, fmt="%2.2g", fontsize=10)
    plt.title(r"$M_A = 200$ GeV")
    plt.xlabel(r"$M_{H^+}$ / GeV")
    plt.ylabel(r"$M_H$ / GeV")
    pdf.savefig()
    plt.close()
