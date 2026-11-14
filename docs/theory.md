# Theory and definitions

## Symmetry-equivalent configurations

A configuration is a set of sites occupied by the dopant (or left vacant). Two configurations are equivalent if a symmetry operation of the crystal, \(\{W|w\}\): \(x \to Wx + w\), maps the sites of one onto the sites of the other. config_enum takes the operations from spglib (for the structure itself, or for the reference structure given with `-ref`), works out how each operation permutes the selected sites, and keeps the first configuration of every set of equivalent ones. The **degeneracy** \(\Omega_m\) of the kept configuration \(m\) is the number of configurations in its set. The degeneracies add up to the total number of combinations, \(\sum_m \Omega_m = \binom{N}{n}\) for \(n\) dopants on \(N\) sites.

Two positions are taken to be the same site if the squared distance of their fractional coordinates is below the tolerance `-tol` (default 10⁻³), which is also the `symprec` given to spglib.

### Sequential doping

In step 2 and later of a sequential doping (`-seq`), the starting structure already contains dopants, and its symmetry is lower. The new dopants are enumerated with the symmetry operations of this doped structure. Two of the resulting structures can still be equivalent through an operation of the original (Step 1) structure that moves the dopants placed earlier — for example, adding B to a structure with A, and adding A to a structure with B. config_enum therefore compares the complete sets of changed sites of the new structures under the symmetry of the Step 1 structure and merges the equivalent ones (their degeneracies are added). The same test prevents equivalent configurations made from different input structures from being picked twice.

With `-ssym ref`, the symmetry operations of the Step 1 structure are used directly in the later steps, as the original script did. This treats sites as equivalent even if they are at different distances from the dopants already present, so it gives fewer configurations but can miss inequivalent ones.

## Configurational thermodynamics

The statistics follow R. Grau-Crespo *et al.*, J. Phys.: Condens. Matter **19**, 256201 (2007). For \(M\) inequivalent configurations with energies \(E_m\) and degeneracies \(\Omega_m\) at temperature \(T\) (`-t`):

- reduced energy: \(\tilde E_m = E_m - T S_m\) with the degeneracy entropy \(S_m = k_B \ln \Omega_m\) (column `E_reduced`);
- probability: \(\tilde P_m = \exp(-\tilde E_m / k_B T) / Z\) with \(Z = \sum_m \exp(-\tilde E_m / k_B T) = \sum_m \Omega_m \exp(-E_m/k_B T)\) (column `Prob.`);
- average energy \(E = \sum_m \tilde P_m E_m\), free energy \(F = -k_B T \ln Z\), and configurational entropy \(S = (E - F)/T\), with the upper limit \(S_\mathrm{max} = k_B \ln \sum_m \Omega_m\) (printed with `-v`; returned by `config_enum.thermo.config_thermo`).

A configuration with a higher degeneracy has a higher configurational entropy and thus a lower reduced energy. The statistics is done separately for each input structure and each number of dopants, since only configurations of the same composition can be compared.

!!! note
    The energies are those of the engine: total energies of the (relaxed) configurations. Vibrational contributions are not included, and the result depends on the cell size, as only configurations within the chosen supercell are sampled.

## Picking configurations (`-stype`)

| `-stype` | Order |
|---|---|
| `E` (default) | lowest reduced energy first |
| `S` | highest space-group number first; lowest reduced energy first among those with the same space group |
| `D` | highest degeneracy first |
| `R` | random |

The first `-spick` configurations are picked. Configurations whose calculation failed are picked only if there are no others.
