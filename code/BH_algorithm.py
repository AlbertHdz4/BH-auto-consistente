""" Alberto Hernandez Lopez """
""" Archivo para implementar el 
    algoritmo de BH - autoconsistente """

from hamiltonian    import *


def BH_algorithm (t, mu, error, tolerance, psi, n_operator, b, b_dagger, identity) : 
    """ Calcula las fluctuaciones y los nuevos estados base
        t           : Parametro t
        mu          : Potencial quimico
        error       : Error inicial
        tolerance   : Tolerancia a alcanzar
        psi         : Ansatz de Psi
        n_operator  : Operador de numero 
        b           : Operador de ascenso
        b_dagger    : Operador de descenso
        identity    : Matriz identidad
        return 
            psi     : Estado final aproximado 
            delta_n : Fluctuaciones
            i       : Iteraciones realizadas
    """
    i = 0
    isAGoodApproximation = True
    
    while (isAGoodApproximation) :
        H           = get_hamiltonian(t, mu, psi, n_operator, \
                        b, b_dagger, identity)
        base_state  = get_base_state(H)
        new_psi     = np.dot(np.transpose(base_state), np.dot(b, base_state))
        error       = np.abs(psi - new_psi)
        psi         = new_psi
        isAGoodApproximation = (error > tolerance)
        i += 1

        # print('error: ', error)
        # print('tolerance: ', tolerance)
        # print('iteration: {0}'.format(i))

    n_squared               = np.dot(n_operator, n_operator)
    base_state_transpose    = np.transpose(base_state)
    squared_n_average       = np.dot(base_state_transpose, np.dot(n_squared, base_state))
    squared_average_n       = (np.dot(base_state_transpose, np.dot(n_operator, base_state))) ** 2
    delta_n                 = squared_n_average - squared_average_n

    return psi, delta_n, i


def parabolic_function (mu_range) :
    """ Generamos la funcion parabolica y = -g^2 + 2xg + g - x^2 - x
        mu_range : Rango del potencial quimico [0, 3]
        return
            parabolic_function : Funcion parabolica
    """
    x = mu_range
    g = np.trunc(mu_range + 1)
    return lambda parabolic_function: -(g) ** 2 + 2 * x * g + g - x ** 2 


def Uzt_BH_algorithm (mu, u, error, tolerance, psi, n_operator, b, b_dagger, identity) : 
    """ Calcula las fluctuaciones y los nuevos estados base
        mu          : Potencial quimico
        error       : Error inicial
        tolerance   : Tolerancia a alcanzar
        psi         : Ansatz de Psi
        n_operator  : Operador de numero 
        b           : Operador de ascenso
        b_dagger    : Operador de descenso
        identity    : Matriz identidad
        return 
            psi     : Estado final aproximado 
            delta_n : Fluctuaciones
            i       : Iteraciones realizadas
    """
    i = 0
    isAGoodApproximation = True
    
    while (isAGoodApproximation) :
        H           = get_reduced_hamiltonian(mu, psi, u, \
                        b, b_dagger, n_operator, identity)
        base_state  = get_base_state(H)
        new_psi     = np.dot(np.transpose(base_state), np.dot(b, base_state))
        error       = np.abs(psi - new_psi)
        psi         = new_psi
        isAGoodApproximation = (error > tolerance)
        i += 1

    n_squared               = np.dot(n_operator, n_operator)
    base_state_transpose    = np.transpose(base_state)
    squared_n_average       = np.dot(base_state_transpose, np.dot(n_squared, base_state))
    squared_average_n       = (np.dot(base_state_transpose, np.dot(n_operator, base_state))) ** 2
    delta_n                 = squared_n_average - squared_average_n

    return psi, delta_n, i



def gen_data_to_plot (t_or_u_range, mu_range, error, tolerance, psi, n_operator, b, b_dagger, identity, steps, option) : 
    """ Generamos los datos para N pasos de la funcion Bose - Hubbard
        t_range     : Rango de t, para nuestro caso va de [0, 1]
        mu_range    : Rango de mu, para nuestro casi va de [0, 3]
        error       : Error inicial
        tolerance   : Tolerancia a alcanzar
        psi         : Ansatz de Psi
        n_operator  : Operador de numero 
        b           : Operador de ascenso
        b_dagger    : Operador de descenso
        identity    : Matriz identidad
        steps       : Pasos para la generacion de los datos
        return 
            psis    : Valores de los estados calculados aproximados
            deltas  : Valores de las fluctuaciones 
    """
    psis    = np.zeros((steps, steps)) # Los estados calculados
    deltas  = np.zeros((steps, steps)) # Fluctuaciones
    iterations = np.zeros((steps, steps))

    if option == 1 : 
        for i in range(steps) :
            for j in range (steps) :
                psis[i, j], deltas[i, j], iterations[i, j] = \
                    BH_algorithm(t_or_u_range[i], mu_range[j], error,\
                        tolerance, psi, n_operator, b, b_dagger, identity)
    else :
        
        for i in range(steps) :
            for j in range (steps) :
                psis[i, j], deltas[i, j], iterations[i, j] = \
                    Uzt_BH_algorithm(mu_range[j], t_or_u_range[i], error,\
                        tolerance, psi, n_operator, b, b_dagger, identity)
                        
    return psis, deltas, iterations