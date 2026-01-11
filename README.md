# LRON_SW

This code provides a Monte Carlo simulation for the O(N) spin model in any spin dimension $N$ and any spatial dimension $d$. More importantly, it can simulate long-range interactions characterized by a power-law decay, expressed as $1/r^{d+\sigma}$.

Here, the Hamiltonian of the $d$-dimensional LR-O(N) model on a hypercubic lattice with periodic boundary conditions (PBCs) is defined as

$$H = -\sum_{i < j} \frac{c(\sigma, L)}{r_{ij}^{d+\sigma}}\mathbf{S}_i \cdot \mathbf{S}_j$$

where $\mathbf{S}_i$ denotes an $N$-component unit vector at site $i$, and $N = L^d$ is the total number of sites. The interactions follow the shortest path on the surface, corresponding to the minimum-image convention. Hence, the summation in the Hamiltonian runs over all spin pairs, resulting in $N(N-1)/2$ interaction terms. The normalization constant $c(\sigma, L)$ is chosen to satisfy

$$\sum_{j=2}^{N}\frac{c(\sigma, L)}{r_{1j}^{d+\sigma}} = 2d$$

where the value $2d$ corresponds to the number of nearest-neighbor (NN) sites. This normalization ensures that as $\sigma \to \infty$, the critical point of the system approaches the standard nearest-neighbor critical temperature $T_{\rm NN}$, while for small $\sigma$, the critical temperature tends to $T_c = 2d/N$.

The update scheme is based on the Clock technique; see the paper: [Phys. Rev. E 99, 010105(R) (2019)](https://doi.org/10.1103/PhysRevE.99.010105).


## How to use

1. Clone the repository:
    ```bash
    git clone https://github.com/Tensofermi/LRON_SW
    cd LRON_SW
    ```

2. Configure model selection in `Makefile`:
    ```bash
    SRC_DIR := src
    # SRC_DIR := src_2D_Ising_v1
    # SRC_DIR := src_2D_Ising_v2
    ```
    Uncomment the desired model and comment out the others. The latter two are designed specifically for the 2D long-range Ising model, but they can also measure the geometric quantity in the FK representation—the wrapping probability—to calculate the critical polynomial $R_p$. The difference between v1 and v2 lies in the method used to measure the wrapping probability: v1 performs a search based on the BFS algorithm, while v2 measures it dynamically during the update process, following the method described in the paper: [Tolson H. Bell et al., 2022, *J. Phys. A: Math. Theor.* 55, 044001](https://iopscience.iop.org/article/10.1088/1751-8121/ac42ab).

3.  Set parameters in `input.txt` and run the simulation:
    ```
    ./run.sh
    ```

4. For advanced simulations, including job submissions on HPC systems, refer to the `/lsub`, `/qsub` and `/data` directories in our related project:  [Zoo_of_Classical_ON_Spin_Model.](https://github.com/Tensofermi/Zoo_of_Classical_ON_Spin_Model)

## Citation
If you use this code or refer to our results in your research, please cite:
- LR-Ising
``` TEXT
@article{xiao2025sak,
  title={On Sak's criterion for statistical models with long-range interaction},
  author={Xiao, Tianning and Liu, Ziyu and Fan, Zhijie and Deng, Youjin},
  journal={arXiv preprint arXiv:2512.04805},
  year={2025}
}
```
- LR-XY
```
@article{xiao2024two,
  title={Two-dimensional xy ferromagnet induced by long-range interaction},
  author={Xiao, Tianning and Yao, Dingyun and Zhang, Chao and Fan, Zhijie and Deng, Youjin},
  journal={Chinese Physics Letters}
}

@article{yghs-p5mg,
  title = {Nonclassical regime of the two-dimensional long-range XY model: A comprehensive Monte Carlo study},
  author = {Yao, Dingyun and Xiao, Tianning and Zhang, Chao and Deng, Youjin and Fan, Zhijie},
  journal = {Phys. Rev. B},
  volume = {112},
  issue = {14},
  pages = {144429},
  numpages = {26},
  year = {2025},
  month = {Oct},
  publisher = {American Physical Society},
  doi = {10.1103/yghs-p5mg},
  url = {https://link.aps.org/doi/10.1103/yghs-p5mg}
}
```

- LR-Heisenberg
```
@article{yao2025spontaneous,
  title={Spontaneous Symmetry Breaking in Two-dimensional Long-range Heisenberg Model},
  author={Yao, Dingyun and Xiao, Tianning and Fan, Zhijie and Deng, Youjin},
  journal={arXiv preprint arXiv:2512.01956},
  year={2025}
}
```
