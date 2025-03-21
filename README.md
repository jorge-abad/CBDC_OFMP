# CBDC and the operational framework of monetary policy

Replication codes for Abad, Nuño and Thomas "[CBDC and the operational framework of monetary policy](https://www.sciencedirect.com/science/article/pii/S0304393225000339)" (Journal of Monetary Economics, 2025)

### Repository Structure

  •	`calibration/` – Contains functions used in the model calibration. \
	•	`functions/` – Contains functions used in the model solution. \
	•	`transitions/` – Contains Dynare files used in the transitional dynamics analysis. \
	•	`graphs/` – Stores output figures. \
	•	`A0_script.m` – This is the main file and should be run to generate the results in the paper.

### Requirements
 
Before running `A0_script.m`, you might need to modify **line 17** to specify the correct path to your Dynare installation.

```matlab
addpath '/Applications/Dynare/6.1-arm64/matlab'
```

### Citation

If you use this code, please cite the paper as follows:
```bibtex
@article{ABAD2025103762,
title = {CBDC and the operational framework of monetary policy},
journal = {Journal of Monetary Economics},
pages = {103762},
year = {2025},
issn = {0304-3932},
doi = {https://doi.org/10.1016/j.jmoneco.2025.103762},
url = {https://www.sciencedirect.com/science/article/pii/S0304393225000339},
author = {Jorge Abad and Galo Nuño and Carlos Thomas},
}
```

### License

This code is released under the **MIT License**. You are free to use, modify, and distribute it, provided that you include proper attribution.  

See the [LICENSE](./LICENSE) file for full details.

### Contact

For questions or issues, please contact [jorgeabad@bde.es](mailto:jorgeabad@bde.es).
