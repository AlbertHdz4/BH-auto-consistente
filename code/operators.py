""" Alberto Hernandez Lopez """
""" Archivo para formar los operadores """

import numpy as np


def get_b_dagger_operator (dimen) :
    """ Forma el operador b daga
        dimen : Dimension del operador b daga
    """
    return np.zeros((dimen, dimen))
    

def get_b_operator (dimen) :
    """ Forma el operador b
        dimen : Dimension del operador b
    """
    return np.zeros((dimen, dimen))


def get_n_operator (dimen) : 
    """ Forma el operador de numero
        dimen : Dimension del operador de numero
    """
    return np.zeros((dimen, dimen))


def gen_b_dagger_b_n_operator (dimen) : 
    """ Genera los operadores b daga, b y n
        dimen : Dimension de los operadores
        return 
            n : Operador de numero 
            b : Operador de descenso
            b_dagger : Operador de ascenso
            identity : Matriz identidad
    """
    b_dagger    = get_b_dagger_operator(dimen)
    b           = get_b_operator(dimen)
    n           = get_n_operator(dimen)
    identity    = np.eye(dimen)
    for i in range(dimen - 1) :
        n[i, i]             = i
        b[i, i + 1]         = np.sqrt(i + 1)
        b_dagger[i + 1, i]  = np.sqrt(i + 1)

    n[dimen - 1, dimen - 1] = dimen - 1
    return (n, b, b_dagger, identity)