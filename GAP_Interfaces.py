from quippy.descriptors import Descriptor
from quippy.potential import Potential
from xml.etree.ElementTree import parse
import os
import numpy as np
import scipy
from ase.data import chemical_symbols
import re

###################################
# PYTHON GAP MODEL IMPLEMENTATION #
###################################

# RegEx to search command line for "Z=, Z1 = or z2= style args"
Z_regex = "(Z|z)[1-9]*\s?=\s?([1-9]+)"
Z_regex = re.compile(Z_regex)


# Main method to calculate cov kernel from cov type flag

def cov_kernel(cov_type, xs, x_cuts, ys, y_cuts, cov_prop):
    '''
    Generate covariance kernel matrix for several kinds of similarity function

    cov_type : int
      Key giving kind of similarity function to use. 1:ard_se, 2:dot_product

    xs: numpy.ndarray
      X values

    '''
    nx = len(x_cuts)
    ny = len(y_cuts)

    if cov_type == 1:  # ard_se
        theta = cov_prop

        K = _ard_se(theta, xs, ys, nx, ny)

    elif cov_type == 2:  # dot_product
        zeta = cov_prop

        K = _dot_product(xs, ys, nx, ny) ** zeta

    for i in range(nx):
        for j in range(ny):
            K[i, j] *= x_cuts[i] * y_cuts[j]
    return K

# Covariance Kernels


def _ard_se(theta, xs, ys, nx, ny):
    K = np.zeros((nx, ny))

    for i in range(nx):
        for j in range(ny):
            K[i, j] = (xs[i, 0] - ys[j, 0])

    K = np.exp(-K**2 / (2 * theta**2))
    return K


def _dot_product(xs, ys, nx, ny):
    K = np.zeros((nx, ny))

    for i in range(nx):
        K[i, :] = xs[i, ...] @ ys.T
    return K


# Python implementation of the GAP model
class GapPy():
    # Initialise GAP model
    def __init__(self, gp_dir, gp_file, sigma):

        self.sigma = sigma

        self.fname = gp_dir + os.sep + gp_file

        # Open XML file
        tree = parse(self.fname)
        self.gp_xml = tree.getroot()

        # GP Label
        self.gp_label = self.gp_xml[0].attrib["label"]

        # Grab descriptor branches
        self.all_desc = self.gp_xml[1][1][:]

        # Isolated atom energies
        isolated_energies = self.gp_xml[1][0][:]
        self.isolated_energies = {}
        for at_type in isolated_energies:
            self.isolated_energies[int(at_type.attrib["Z"])] = float(
                at_type.attrib["value"])

        self.num_desc = len(self.all_desc)

        # Space to dump descriptor params
        self.descs = []
        self.desc_types = []
        self.desc_cmds = []
        self.desc_names = []
        self.nsparses = []
        self.alphas = []
        self.cuts = []
        self.cov_types = []
        self.cov_props = []
        self.deltas = []

        # Grab details for each descriptor
        for desc in self.all_desc:
            descriptor = [child for child in desc if "descriptor" in child.tag]

            alphas = np.array([float(child.attrib["alpha"])
                              for child in desc if "sparseX" in child.tag])
            cuts = np.array([float(child.attrib["sparseCutoff"])
                            for child in desc if "sparseX" in child.tag])
            delta = float(desc.attrib["signal_variance"])
            cov_type = int(desc.attrib["covariance_type"])

            cmd = descriptor[0].text

            if cov_type == 1:  # ard_se
                cov_prop = float(
                    [child.text for child in desc if "theta" in child.tag][0])
            elif cov_type == 2:  # dot_product
                cov_prop = float(desc.attrib["zeta"])
            else:
                print(f"Cov type {cov_type} not understood")
                cov_prop = 0

            desc_type = cmd.split(" ")[0]

            # Find chemical species the descriptor acts on
            Zs = Z_regex.findall(cmd)
            desc_Zs = [Z[1] for Z in Zs]

            # Generate a helpful name for the descriptor
            desc_name = "".join([chemical_symbols[int(Z)]
                                for Z in desc_Zs]) + " " + desc_type

            self.descs.append(Descriptor(descriptor[0].text))
            self.desc_cmds.append(cmd)
            self.desc_types.append(desc_type)
            self.desc_names.append(desc_name)
            self.nsparses.append(len(alphas))
            self.alphas.append(alphas)
            self.cuts.append(cuts)
            self.cov_types.append(cov_type)
            self.cov_props.append(cov_prop)
            self.deltas.append(delta)

        # Load sparse points for each descriptor from .sparseX files
        self.sparse_files = sorted([file for file in os.listdir(
            gp_dir) if gp_file in file and "sparseX" in file and self.gp_label in file])

        self.sparseXs = [np.loadtxt(gp_dir + os.sep + file).reshape(
            (self.nsparses[i], -1)) for i, file in enumerate(self.sparse_files)]

        # Generate sparseX-sparseX covariance Kernels (K_xx) & Cholesky decomps
        self.Ks = []
        self.Ls = []

        for i in range(self.num_desc):
            cov_type = self.cov_types[i]
            cov_prop = self.cov_props[i]

            delta = self.deltas[i]

            sparseX = self.sparseXs[i]
            cuts = self.cuts[i]

            K = cov_kernel(cov_type, sparseX, cuts, sparseX,
                           cuts, cov_prop) * delta**2
            K += self.sigma**2 * np.eye(K.shape[0])

            self.Ks.append(K)

            L = np.linalg.cholesky(K)
            self.Ls.append(L)

    def predict_energy(self, atoms):
        '''
        Predict total energy posterior mean and variance
        '''

        nats = len(atoms)

        # Energy per atom and covariance matrix of atomic energies
        E_tot = np.zeros(nats)
        E_cov = np.zeros((nats, nats))

        # Account for isolated atom energies
        E_tot = np.array([self.isolated_energies[Z]
                         for Z in atoms.get_atomic_numbers()])

        # Loop over each descriptor
        for i, desc in enumerate(self.descs):
            # Get kind of cov kernel (dot_product, ard_se, ...)
            cov_type = self.cov_types[i]
            # Get kernel hyperparameters (theta, xi, ...)
            cov_prop = self.cov_props[i]

            # Energy scaling factor for descriptor
            delta = self.deltas[i]

            # Weights
            alphas = self.alphas[i]

            # Apply descriptor to model
            result = desc.calc(atoms, grad=True)

            # Descriptor vectors found in structure
            x_star = result["data"]

            x_star_cut = result["covariance_cutoff"]

            # Atomic indices
            indices = result["ci"].reshape((len(x_star), -1))
            npasses = indices.shape[-1]

            indices -= 1  # Convert from Fortran indexing!

            x = self.sparseXs[i]
            x_cut = self.cuts[i]

            # Generate all cov matrices
            L = self.Ls[i]  # Cholesky decomp of (K_x_x + Sigma)

            # Construct all useful matrices
            K_xstar_x = cov_kernel(
                cov_type, x_star, x_star_cut, x, x_cut, cov_prop) * delta**2

            K_x_xstar = K_xstar_x.T

            K_xstar_xstar = cov_kernel(
                cov_type, x_star, x_star_cut, x_star, x_star_cut, cov_prop) * delta**2

            # Energy prediction
            # sol computes [K_x_x + Sigma]^-1 via Cholesky decomp
            sol = scipy.linalg.solve_triangular(L, K_x_xstar, lower=True).T

            # Expressions for the descriptor energies and covariance matrix
            E_desc = K_xstar_x @ alphas
            cov_desc = K_xstar_xstar - sol @ sol.T

            for j in range(npasses):
                for i in range(indices.shape[0]):
                    # Distribute descriptor energies onto atoms
                    # Multiple passes required, as N-body potentials
                    # Have one descriptor energy split between N atoms
                    # EG: 2b bond gives 1/2 bond energy to each atom
                    E_tot[indices[i, j]] += E_desc[i] / npasses

                    for k in range(indices.shape[0]):
                        E_cov[indices[i, j], indices[k, j]
                              ] += cov_desc[i, k] / npasses

        return E_tot, E_cov


############################################
# INTERFACE FOR QUIP IMPLEMENTATION OF GAP #
############################################

def GapQUIP():
    fpath = os.path.dirname(os.path.abspath(__file__)) + os.sep
    gap = Potential(param_filename=fpath +
                    f"GAP/InP_GAP.xml", calc_args='local_gap_variance')
    return gap
