# Open-channel hydraulics: normal and critical depth

Jupyter notebooks that compute the normal depth (Manning's equation) and critical depth
(Froude number = 1) of prismatic channel sections with the Newton–Raphson method, with an
interactive cross-section plot.

## Notebooks

| Notebook | Section | Colab | Binder |
| --- | --- | --- | --- |
| [Ch2Depth_Rectangle.ipynb](normal_critical_depth/Ch2Depth_Rectangle.ipynb) | Rectangle | <a href="https://colab.research.google.com/github/AaronOET/ntouoet_openchannel_1151/blob/main/normal_critical_depth/Ch2Depth_Rectangle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a> | <a href="https://mybinder.org/v2/gh/AaronOET/ntouoet_openchannel_1151/main?labpath=normal_critical_depth/Ch2Depth_Rectangle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a> |
| [Ch2Depth_Trapezoid.ipynb](normal_critical_depth/Ch2Depth_Trapezoid.ipynb) | Trapezoid | <a href="https://colab.research.google.com/github/AaronOET/ntouoet_openchannel_1151/blob/main/normal_critical_depth/Ch2Depth_Trapezoid.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a> | <a href="https://mybinder.org/v2/gh/AaronOET/ntouoet_openchannel_1151/main?labpath=normal_critical_depth/Ch2Depth_Trapezoid.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a> |
| [Ch2Depth_Triangle.ipynb](normal_critical_depth/Ch2Depth_Triangle.ipynb) | Triangle | <a href="https://colab.research.google.com/github/AaronOET/ntouoet_openchannel_1151/blob/main/normal_critical_depth/Ch2Depth_Triangle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a> | <a href="https://mybinder.org/v2/gh/AaronOET/ntouoet_openchannel_1151/main?labpath=normal_critical_depth/Ch2Depth_Triangle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a> |
| [Ch2Depth_Circle.ipynb](normal_critical_depth/Ch2Depth_Circle.ipynb) | Circle (part-full pipe) | <a href="https://colab.research.google.com/github/AaronOET/ntouoet_openchannel_1151/blob/main/normal_critical_depth/Ch2Depth_Circle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a> | <a href="https://mybinder.org/v2/gh/AaronOET/ntouoet_openchannel_1151/main?labpath=normal_critical_depth/Ch2Depth_Circle.ipynb" target="_blank" rel="noopener noreferrer"><img src="https://mybinder.org/badge_logo.svg" alt="Binder"></a> |

Each notebook lists the geometric functions of its section (Chow, 1959, Table 2.1), checks
the derivative term against a finite difference, works an example and ends with sliders
for the geometry, bed slope $S_0$, Manning's $n$ and discharge $Q$.

## Modules

- [channel_sections.py](normal_critical_depth/channel_sections.py): section geometry and the normal- and critical-depth solvers
- [channel_plot.py](normal_critical_depth/channel_plot.py): cross-section plot and ipywidgets controls

## Running locally

```bash
pip install -r requirements.txt
jupyter lab
```

## References

- Chow, V. T., Maidment, D. R., & Mays, L. W. (1988). *Applied Hydrology*. New York: McGraw-Hill.
- Kuiry, S. N., & Sen, D. (2020). *Modelling Hydrology, Hydraulics and Contaminant Transport
  Systems in Python*. Boca Raton, FL: CRC Press (Taylor & Francis Group).
