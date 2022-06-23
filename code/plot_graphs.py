""" Alberto Hernandez Lopez """
""" Archivo para graficar los datos """

import matplotlib.pyplot as plt


def graph_ztU_vs_muU (funct_data, mu_range) :
    plt.plot(funct_data, mu_range, "b")
    plt.xlabel('zt/U', size = 12)
    plt.ylabel(r'$\mu/U$', size = 12)
    plt.grid()
    plt.legend()
    plt.show()


def graph_warm_psi_squared_ztU_vs_muU (funct_data, mu_range, t_values, transpose_psis) : 
    plt.plot(funct_data, mu_range, "b")
    plt.imshow(transpose_psis, extent = t_values, aspect = 'auto', cmap = 'magma', origin = 'lower')
    plt.xlabel('zt/U',size = 12)
    plt.ylabel(r'$\mu/U$',size = 12)
    plt.xlim([0, 1.])
    plt.colorbar().set_label('|psi|^2', size = 12)
    plt.legend()
    plt.show()


def graph_warm_iterations_ztU_vs_muU (iterations, parameters) :
    plt.imshow(iterations, extent = parameters, aspect = 'auto', cmap = 'plasma', origin = 'lower')
    plt.xlabel('zt/U',size = 12)
    plt.ylabel(r'$\mu/U$',size = 12)
    plt.xlim([0, 1.])
    plt.colorbar().set_label('iterations', size = 12)
    plt.show()


def graph_warm_fluctuations_ztU_vs_muU (fluctuations, parameters) : 
    plt.imshow(fluctuations, extent = parameters, aspect = 'auto', cmap = 'plasma', origin = 'lower')
    plt.xlabel('zt / U',size = 12)
    plt.ylabel(r'$\mu/U$',size = 12)
    plt.xlim([0, 1.])
    plt.colorbar().set_label('Fluctuations', size = 12)
    plt.show()


def graph_g1_g2_g3 (u_range, f1_g1_data, f2_g1_data, f1_g2_data, f2_g2_data, f1_g3_data, f2_g3_data) :
    plt.plot(u_range, f1_g1_data, "g")
    plt.plot(u_range, f2_g1_data, "g")

    plt.plot(u_range, f1_g2_data, "g")
    plt.plot(u_range, f2_g2_data, "g")

    plt.plot(u_range, f1_g3_data, "g")
    plt.plot(u_range, f2_g3_data, "g")

    plt.ylabel(r'$\mu/zt$', size = 12)
    plt.xlabel('U/zt', size = 12)
    plt.grid()
    plt.legend()
    plt.show()


def graph_warm_psis_g1_g2_g3 (psis, parameters, u_range, f1_g1_data, f2_g1_data, f1_g2_data, f2_g2_data, f1_g3_data, f2_g3_data) :
    plt.imshow(psis, extent = parameters, aspect = 'auto', cmap = 'inferno', origin = 'lower')
    plt.plot(u_range, f1_g1_data, "g")
    plt.plot(u_range, f2_g1_data, "g")

    plt.plot(u_range, f1_g2_data, "g")
    plt.plot(u_range, f2_g2_data, "g")

    plt.plot(u_range, f1_g3_data, "g")
    plt.plot(u_range, f2_g3_data, "g")

    plt.ylabel(r'$\mu/zt$', size = 12)
    plt.xlabel('U/zt', size = 12)
    plt.colorbar().set_label(r'$|\psi|^{2}$',size = 12)
    plt.legend()
    plt.show()


def grap_warm_iterations_g1_g2_g3 (iterations, parameters) : 
    plt.imshow(iterations, extent = parameters, \
                aspect = 'auto', cmap = 'plasma', \
                    origin = 'lower')
    plt.ylabel(r'$\mu/zt$', size = 12)
    plt.xlabel('U/zt', size = 12)
    plt.colorbar().set_label('iteraciones', size = 12)
    plt.show()


def grap_warm_fluctuations_g1_g2_g3 (fluctuations, parameters) : 
    plt.imshow(fluctuations, extent = parameters, aspect = 'auto', cmap='plasma', origin='lower')
    plt.ylabel(r'$\mu/zt$', size = 12)
    plt.xlabel('U/zt', size = 12)
    plt.colorbar().set_label('fluctuaciones', size = 12)
    plt.show()