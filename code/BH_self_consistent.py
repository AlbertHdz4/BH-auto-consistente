""" Alberto Hernandez Lopez """
""" Archivo principal para ejecutar 
    el BH autoconsistente """

from hamiltonian    import *
from operators      import *
from plot_graphs    import *
from BH_algorithm   import *


if __name__ == '__main__' :
    dimension   = 30
    tolerance   = 0.0001
    error       = 1.
    psi         = 1. # Ansatz Psi0
    mu          = 1. # potencial quimico
    t           = 1. # Amplitud de tunelamiento

    (n_operator, b, b_dagger, identity) = \
        gen_b_dagger_b_n_operator(dimension)

    # Test - Para solo un caso
    # psi, delta_n, iterations = BH(t, mu, error, tolerance, psi, n_operator, b, b_dagger, identity)

    N = 350
    t_range     = np.linspace(0, 1, N)
    mu_range    = np.linspace(0, 3, N)
    psis, deltas_n, iterations = gen_data_to_plot(t_range, mu_range, error, tolerance, psi, \
                        n_operator, b, b_dagger, identity, N, option = 1)

    print('psis: ', psis)
    print('deltas_n: ', deltas_n)
    
    x = mu_range
    g = np.trunc(mu_range + 1) # Parte entero, funcion piso

    # Funcion parabolica : y = -g^2 + 2xg + g - x^2 - x
    parabolic_funct = lambda x : (-(g) ** 2 + 2 * x * g + g - x ** 2 - x) / (1 + x)

    funct_data = parabolic_funct(mu_range)
    graph_ztU_vs_muU(funct_data, mu_range)

    parameters = [t_range[0], t_range[N - 1], mu_range[0], mu_range[N - 1]]

    transpose_psis = np.transpose(psis)
    graph_warm_psi_squared_ztU_vs_muU(funct_data, mu_range, \
                                parameters, transpose_psis)

    iterations_transpose = np.transpose(iterations)
    graph_warm_iterations_ztU_vs_muU(iterations_transpose, parameters)

    transpose_fluctuations = np.transpose(deltas_n)
    graph_warm_fluctuations_ztU_vs_muU(transpose_fluctuations, parameters)


    print('******************** Para el caso U/zt ********************') 
    u           = 1. # Fuerza de interacción
    error       = 1.
    psi         = 1. # Ansatz Psi0
    mu          = 1. # potencial quimico
    
    (n_operator, b, b_dagger, identity) = \
        gen_b_dagger_b_n_operator(dimension)

    psi, delta_n, iterations = \
        Uzt_BH_algorithm(mu, u, error, tolerance, psi, \
                        n_operator, b, b_dagger, identity)
    print('psi: ', psi)
    print('delta_n: ', delta_n)
    print('iterations: ', iterations)

    u_range  = np.linspace(0, 20, N)
    mu_range = np.linspace(0, 50, N)
    psis, deltas_n, iterations = gen_data_to_plot(u_range, mu_range, error, tolerance, psi, \
                        n_operator, b, b_dagger, identity, N, option = 2)

    f1 = lambda u, g: 0.5 * (u * (2 * g - 1) - 1) + \
                        0.5 * np.sqrt(u ** 2 - (2 * u * (2 * g + 1)) + 1)

    f2 = lambda u, g: 0.5 * (u * (2 * g - 1) - 1) - \
                        0.5 * np.sqrt(u ** 2 - (2 * u * (2 * g + 1)) + 1)

    graph_g1_g2_g3(u_range, f1(u_range, 1), f2(u_range, 1), \
                    f1(u_range, 2), f2(u_range, 2), \
                    f1(u_range, 3), f2(u_range, 3))


    parameters = [u_range[0], u_range[N - 1], mu_range[0], mu_range[N - 1]]
    graph_warm_psis_g1_g2_g3(np.transpose(psis), parameters, \
                    u_range, f1(u_range, 1), f2(u_range, 1), \
                    f1(u_range, 2), f2(u_range, 2), \
                    f1(u_range, 3), f2(u_range, 3))
    
    grap_warm_iterations_g1_g2_g3(np.transpose(iterations), parameters)

    grap_warm_fluctuations_g1_g2_g3(np.transpose(deltas_n), parameters)
  