"""Calculate factors for hyperfine structure in single-electron atoms."""

import json
import warnings
from itertools import product

import numpy as np

from clerq.angmom import calc_clebsch_hf

# Zeeman conventions. The shift of a sublevel in a DC field B along the
# quantisation axis is  E = +g_F * m_F * mu_B * B  (the convention of Steck's
# alkali data), where g_F is dimensionless and positive for an electron-like
# moment: the stretched state m_F = +F of Rb87 5S1/2 F=2 (g_F = +1/2) moves UP
# in energy with B > 0. This is the sign used in Steck's tables; some other
# references (e.g. Budker et al.) define the magnetic moment as -g_F mu_B m_F,
# giving g_F of the opposite sign, so check the convention when importing values.
MU_B_OVER_H_MHZ_PER_GAUSS = 1.399624493
"""Bohr magneton over Planck's constant, in MHz per gauss (CODATA 2018)."""

G_S = 2.00231930436
"""Electron spin g-factor (CODATA 2018)."""


def gauss_to_natural(b_gauss: float) -> float:
    """Convert a field in gauss to natural units, ``mu_B * B / h`` in MHz.

    Natural units are those in which a sublevel's shift is
    ``g_F * m_F * B`` in the same (cyclic) frequency unit as the system's
    ``energies``. Use this helper when the energies are in MHz (time in
    microseconds).
    """

    return b_gauss * MU_B_OVER_H_MHZ_PER_GAUSS


def lande_g_J(L: float, J: float, S: float = 0.5, g_L: float = 1.0, g_S: float = G_S):
    """Lande g-factor g_J of a fine-structure level, in units of mu_B."""

    if J == 0:
        return 0.0
    ls = J * (J + 1)
    return g_L * (ls - S * (S + 1) + L * (L + 1)) / (2 * ls) + g_S * (
        ls + S * (S + 1) - L * (L + 1)
    ) / (2 * ls)


def lande_g_F(I: float, J: float, F: float, g_J: float, g_I: float = 0.0):  # noqa: E741
    """Hyperfine g-factor g_F of an F level, in units of mu_B.

    ``g_F = g_J [F(F+1) - I(I+1) + J(J+1)] / (2F(F+1))
    + g_I [F(F+1) + I(I+1) - J(J+1)] / (2F(F+1))``, with ``g_I`` the (small)
    nuclear g-factor already expressed in units of mu_B, as tabulated by Steck
    (e.g. -0.000995 for Rb87).
    """

    if F == 0:
        return 0.0
    ff = F * (F + 1)
    return (
        g_J * (ff - I * (I + 1) + J * (J + 1)) + g_I * (ff + I * (I + 1) - J * (J + 1))
    ) / (2 * ff)


class _JsonMixin:
    """Mixin that provides to_json_str() for classes with get_json_dict()."""

    def to_json_str(self) -> str:
        """Return a JSON string representation of this object."""
        return json.dumps(
            self.get_json_dict(), indent=2, separators=None, sort_keys=True
        )


class Atom1e(_JsonMixin):
    """Represents an atom object and contains a list of J levels."""

    def __init__(self, element=None, isotope=None):

        self.element = element
        self.isotope = isotope
        self.F_levels = []

    def __repr__(self):

        return self.to_json_str()

    def add_F_level(self, F_level):

        self.F_levels.append(F_level)

    def get_F_level_idx_map(self):
        """Maps each element of the mF_list to a F_level index."""

        F_level_idx_map = []
        for i, F_level in enumerate(self.F_levels):
            for _ in F_level.mF_levels:
                F_level_idx_map.append(i)
        return F_level_idx_map

    def get_mF_list(self):
        """Unnests the mF levels into a single list of dicts."""

        mF_list = []
        for F_level in self.F_levels:
            item_I = F_level.I
            item_J = F_level.J
            item_F = F_level.F
            for mF_level in F_level.mF_levels:
                mF_list.append(
                    {
                        "I": item_I,
                        "J": item_J,
                        "F": item_F,
                        "mF": mF_level.mF,
                        "energy": mF_level.energy,
                    }
                )
        return mF_list

    def get_g_F_list(self) -> list[float]:
        """Returns the g_F factor of each mF level, in units of mu_B.

        Raises:
            ValueError: if an F level has no g_F (give its ``g_F``, ``g_J`` or
                ``L`` when constructing the ``LevelF``).
        """

        g_F_list = []
        for i, F_level in enumerate(self.F_levels):
            if F_level.g_F is None:
                raise ValueError(
                    f"F level {i} (J={F_level.J}, F={F_level.F}) has no g_F. Give "
                    "g_F, g_J or L to its LevelF."
                )
            g_F_list.extend([F_level.g_F] * len(F_level.mF_levels))
        return g_F_list

    def get_zeeman_levels(self) -> dict[str, list[float]]:
        """Returns the per-level ``g_F`` and ``m_F`` lists for a ``b_field``.

        Merge into the system's ``b_field`` block, e.g.
        ``{"magnitude": 1.0, **atom.get_zeeman_levels()}``.
        """

        return {
            "g_F": self.get_g_F_list(),
            "m_F": [float(m["mF"]) for m in self.get_mF_list()],
        }

    def get_num_mF_levels(self):
        """Returns the total number of mF levels in all J levels."""

        return len(self.get_mF_list())

    def get_energies(self):
        """Returns a list of energies for each mF level in all J levels."""

        return [mF_level["energy"] for mF_level in self.get_mF_list()]

    def get_coupled_levels(
        self, F_level_idxs_a: list[int], F_level_idxs_b: list[int]
    ) -> list[list[int]]:
        """Return all mF-level index pairs between two sets of F levels.

        Args:
            F_level_idxs_a: Indices into self.F_levels for the first group.
            F_level_idxs_b: Indices into self.F_levels for the second group.

        Returns:
            List of [i, j] pairs where i is an mF index from group a and
            j is an mF index from group b.

        Note:
            Selection rules (e.g. forbidding J=J' transitions) are not
            enforced here. Callers are responsible for passing physically
            valid level index sets.
        """

        F_level_idx_map = self.get_F_level_idx_map()
        a_levels = [i for i, idx in enumerate(F_level_idx_map) if idx in F_level_idxs_a]
        b_levels = [i for i, idx in enumerate(F_level_idx_map) if idx in F_level_idxs_b]
        # product returns iterator of tuples, convert to list of lists
        return [list(i) for i in product(a_levels, b_levels)]

    def get_clebsch_hf_factors(
        self, F_level_idxs_a: list[int], F_level_idxs_b: list[int], q: int
    ) -> np.ndarray:
        """Returns the Clebsch-Gordan coefficients for the hyperfine
        transition dipole matrix elements of each coupled level pair.

        Args:
            F_level_idxs_a: Indices of the F levels the transition is from
                (lower level).
            F_level_idxs_b: Indices of the F levels the transition is to
                (upper level).
            q: The field polarisation. Choose from [-1, 0, 1].

        Returns:
            Array of factors, one for each pair from ``get_coupled_levels``.

        """

        mF_list = self.get_mF_list()
        coupled_levels = self.get_coupled_levels(F_level_idxs_a, F_level_idxs_b)
        factors = np.empty(len(coupled_levels))
        for i, cl in enumerate(coupled_levels):
            a = mF_list[cl[0]]
            b = mF_list[cl[1]]
            factors[i] = calc_clebsch_hf(
                J_a=a["J"],
                I_a=a["I"],
                F_a=a["F"],
                mF_a=a["mF"],
                J_b=b["J"],
                I_b=b["I"],
                F_b=b["F"],
                mF_b=b["mF"],
                q=q,
            )
        return factors

    def get_clebsch_hf_factors_iso(
        self, F_level_idxs_a: list[int], F_level_idxs_b: list[int]
    ) -> np.ndarray:
        """Returns the Clebsch-Gordan coefficients for the hyperfine
        transition dipole matrix elements of each coupled level pair, for
        an isotropic field.

        Args:
            F_level_idxs_a: Indices of the F levels the transition is from
                (lower level).
            F_level_idxs_b: Indices of the F levels the transition is to
                (upper level).

        Returns:
            Array of factors, one for each pair from ``get_coupled_levels``.

        Notes:
            - An isotropic field is a field with equal components in all three
              possible polarizations.
            - Any given polarisation of the field only interacts with one of the
              three components of the dipole moment, so it is appropriate to
              average over the couplings (i.e. factor 1/3) rather than sum.
        """

        return self.get_decay_factors(F_level_idxs_a, F_level_idxs_b) / np.sqrt(3.0)

    def get_decay_factors(
        self, F_level_idxs_a: list[int], F_level_idxs_b: list[int]
    ) -> np.ndarray:
        """Returns the factors for the collapse operators of each hyperfine
        coupled level pair.

        Args:
            F_level_idxs_a: Indices of the F levels the transition is from
                (lower level).
            F_level_idxs_b: Indices of the F levels the transition is to
                (upper level).

        Returns:
            Array of factors, one for each pair from ``get_coupled_levels``.

        Notes:
            This is equivalent to the clebsch_hf_factors for all
            polarisations as decay photons are of all polarisations. Note
            that for any coupled levels pair there will be only one
            non-zero factor to sum.
        """

        return (
            self.get_clebsch_hf_factors(F_level_idxs_a, F_level_idxs_b, q=-1)
            + self.get_clebsch_hf_factors(F_level_idxs_a, F_level_idxs_b, q=0)
            + self.get_clebsch_hf_factors(F_level_idxs_a, F_level_idxs_b, q=1)
        )

    def get_strength_factor(
        self, F_level_idx_lower, F_level_idx_upper, mF_level_lower_idx=0
    ):
        """Relative hyperfine transition strength factors.
            Equal for each ground state (mF_level_lower_idx), so this parameter
            never needs to be set, just used for testing that claim.

        Notes:
            - Sum of the matrix elements from a single ground-state sublevel
              to the levels in a particular F' energy level.
            - The sum S_{FF'} is independent of the ground state sublevel chosen.
            - The sum of S_{FF'} over upper F levels should be 1.

        Refs:
            [0]: https://steck.us/alkalidata/rubidium87numbers.pdf

        """

        facts_qm1 = self.get_clebsch_hf_factors(
            [F_level_idx_lower], [F_level_idx_upper], q=-1
        )
        facts_q0 = self.get_clebsch_hf_factors(
            [F_level_idx_lower], [F_level_idx_upper], q=0
        )
        facts_qp1 = self.get_clebsch_hf_factors(
            [F_level_idx_lower], [F_level_idx_upper], q=1
        )

        cl = self.get_coupled_levels([F_level_idx_lower], [F_level_idx_upper])

        idx_map = self.get_F_level_idx_map()
        lower_mF_levels = [
            i for i, idx in enumerate(idx_map) if idx == F_level_idx_lower
        ]
        lower_mF_level = lower_mF_levels[mF_level_lower_idx]
        coupled = [lower_mF_level in i for i in cl]

        factor_sq_sum = (
            np.sum(facts_qm1[coupled] ** 2)
            + np.sum(facts_q0[coupled] ** 2)
            + np.sum(facts_qp1[coupled] ** 2)
        )

        return factor_sq_sum

    def get_json_dict(self):

        json_dict = {
            "element": self.element,
            "isotope": self.isotope,
            "F_levels": [i.get_json_dict() for i in self.F_levels],
        }
        return json_dict


class LevelF(_JsonMixin):
    """Represents an F hyperfine structure level and holds its magnetic
        sublevels mF.

    Attributes:
        F (float): Total atomic angular momentum number F.
        energy (float): Energy of the level.
        mf_levels (list, length 2F+1): Energies of 2F+1 hyperfine
            sublevels.

    Notes:
        - The magnitude of F can take values in the range
            `|J - I| <= F <= J + I`. A ``ValueError`` is raised if this is
            not satisfied.
        - The mF_energies are set _relative_ to the F level energy.
    """

    def __init__(
        self,
        I,  # noqa: E741
        J,
        F,
        energy=0.0,
        mF_energies=None,
        g_F=None,
        g_J=None,
        g_I=None,
        L=None,
    ):
        """
        Args:
            F (float): Total atomic angular momentum number F.
            energy (float): Energy of the level.
            mf_energies (list(LevelMF), length 2F+1): List of 2F+1 hyperfine
            sublevels.
            g_F (float): The level's g-factor in units of mu_B. If not given
                it is calculated from ``g_J`` (or from ``L``) and ``g_I``;
                if none of these is given it stays ``None``, and the level
                cannot be used with a B-field.
            g_J (float): Fine-structure Lande factor. Calculated from ``L``
                (with S = 1/2) if not given.
            g_I (float): Nuclear g-factor in units of mu_B (Steck's
                convention, e.g. -0.000995 for Rb87). If ``g_F`` is calculated
                and ``g_I`` is not given, the nuclear term is neglected (0)
                and a warning is issued; pass ``g_I=0.0`` to silence it.
            L (float): Orbital angular momentum, for calculating ``g_J``.
        """

        self.I = I
        self.J = J
        self.F = F
        self.energy = energy

        F_min = abs(J - I)
        F_max = J + I
        if not (F_min <= F <= F_max):
            raise ValueError(
                f"F={F} is outside the allowed range |J-I| <= F <= J+I "
                f"(|{J}-{I}| = {F_min:.1f}, {J}+{I} = {F_max:.1f})."
            )

        if g_F is None:
            if g_J is None and L is not None:
                g_J = lande_g_J(L=L, J=J)
            if g_J is not None:
                if g_I is None:
                    warnings.warn(
                        "g_I not given: the nuclear term of g_F is neglected "
                        "(g_I=0), which shifts g_F by about 1e-3. Pass g_I "
                        "(e.g. -0.0009951 for Rb87) or g_I=0.0 to silence this.",
                        UserWarning,
                        stacklevel=2,
                    )
                    g_I = 0.0
                g_F = lande_g_F(I=I, J=J, F=F, g_J=g_J, g_I=g_I)
        self.g_F = g_F

        self.build_mF_levels(mF_energies)

    def __repr__(self):
        return self.to_json_str()

    def build_mF_levels(self, mF_energies=None):
        """Builds the mF sublevels of the F level.

        Args:
            mF_energies  (list, length 2F+1): Energies of the 2F+1 hyperfine
                sublevels.
        Returns:
            self.mF_levels
        """
        self.mF_levels = []
        mF_range = self.get_mF_range()
        if not mF_energies:
            mF_energies = [0.0 for i in mF_range]
        if len(mF_energies) != len(mF_range):
            raise ValueError("mF_energies is not the correct length.")
        for i, mF in enumerate(mF_range):
            self.mF_levels.append(LevelMF(mF=mF, energy=self.energy + mF_energies[i]))

        return self.mF_levels

    def get_mF_range(self):
        """Returns a range representing the m_F angular momentum sublevels
        [-F, -F+1, ..., F-1, F].
        """

        return np.arange(-self.F, self.F + 1, dtype=float)

    def get_json_dict(self):

        json_dict = {
            "I": self.I,
            "J": self.J,
            "F": self.F,
            "energy": self.energy,
            "mF_levels": [mF.__dict__ for mF in self.mF_levels],
        }
        if self.g_F is not None:
            json_dict["g_F"] = self.g_F
        return json_dict


class LevelMF(_JsonMixin):
    """Represents an m_F hyperfine sublevel.

    Attributes:
        mF (float): Angular momentum number m_F
        energy (float): The energy of the level
    """

    def __init__(self, mF, energy=0.0):
        self.mF = mF
        self.energy = energy

    def __repr__(self):
        return "<LevelMF :: %s>" % self.__dict__

    def get_json_dict(self):
        return self.__dict__
