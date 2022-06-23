""" Alberto Hernandez Lopez """
""" Archivo para calcular el hamiltoniano """

from scipy.linalg   import eigh
from operators      import *

def get_hamiltonian (t, mu, psi, n_operator, b, b_dagger, identity) : 
    """ Obtiene el hamiltoniano del sistema
        mu          : Potencial quimico
        psi         : Estado del sistema
        n_operator  : Operador de numero 
        b           : Operador de descenso 
        b_dagger    : Operador de ascenso
        identity    : Matriz identidad    
        return 
            H : Hamiltoniano del sistema
    """
    kinetic_term    = -t * psi * (b_dagger + b) + t * (psi ** 2) * identity
    H               = kinetic_term + 0.5 * np.dot(n_operator, (n_operator - identity)) - \
                    mu * n_operator
    return H


def get_reduced_hamiltonian (mu, psi, u, b, b_dagger, n_operator, identity) :
    """ Obtiene el Hamiltoniano reducido de un solo sitio
        mu          : Potencial quimico
        psi         : Estado del sistema
        n_operator  : Operador de numero 
        b           : Operador de descenso 
        b_dagger    : Operador de ascenso
        identity    : Matriz identidad    
        return 
            H : Hamiltoniano reducido 
    """
    H = - psi * (b_dagger + b) + (psi ** 2) * identity + \
        0.5 * u * np.dot(n_operator, (n_operator - identity)) - \
        mu * n_operator
    return H



def get_base_state (hamiltonian) :
    """ Obtiene el estado base
        hamiltonian : Hamiltoniano del cual se obtendra el estado base
        return 
            base_state : Regresa el estado base del sistema
    """
    eigenvalues, eigenvectors = eigh(hamiltonian)
    base_state = eigenvectors[:, 0]
    return base_state