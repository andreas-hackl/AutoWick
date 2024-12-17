import numpy as np

# Gamma matrix definition of GPT https://github.com/lehner/gpt

gamma = {
    0: np.array( [[0, 0, 0, 1j], [0, 0, 1j, 0], [0, -1j, 0, 0], [-1j, 0, 0, 0]], dtype=np.cdouble),
    1: np.array( [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]], dtype=np.cdouble),
    2: np.array( [[0, 0, 1j, 0], [0, 0, 0, -1j], [-1j, 0, 0, 0], [0, 1j, 0, 0]], dtype=np.cdouble),
    3: np.array( [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.cdouble),
    5: np.diagflat([1, 1, -1, -1]).astype(dtype=np.cdouble),
    'I': np.diagflat([1, 1, 1, 1]).astype(dtype=np.cdouble)
}

parity = {
    'Pp': 1/2 * (gamma['I'] + gamma[3]),
    'Pn': 1/2 * (gamma['I'] - gamma[3]),
}

C = 1j * gamma[1] @ gamma[3]

nucleon = {
    'Cg5': C @ gamma[5],
    'ig3Cg5': 1j * gamma[3] @ C @ gamma[5],
}

pion = {
    'g5': gamma[5],
    'g3g5': gamma[3] @ gamma[5]
}
