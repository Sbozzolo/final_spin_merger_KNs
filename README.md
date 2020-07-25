[![GPLv3
license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

# final_spin_merger_KNs

Estimate the spin of the final black hole of the merger of two Kerr-Newman black
holes. Implements the method described in arxiv:1706.06519.

# Installation

Clone this repository:
```bash
   git clone https://github.com/Sbozzolo/two_charged_body_problem.git
```
Move into the folder and install with pip:

```bash
   cd two_charged_body_problem.git && pip3 install --user
```

# Usage

Just run `estimate_spins.py`. The available options are:
```
--m1: First mass
--m2: Second mass
--M: Total mass
--q: Mass ratio
--tolerance: Error in determining the spin
--lambda1: First charge/mass ratio
--lambda2: Second charge/mass ratio
```



