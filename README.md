# [Overlapping finite element meshes in AMORE](https://doi.org/10.1016/j.advengsoft.2020.102791)

## Usage

```
python3 femMain.py
SimpleOFE3
```

## Numerical Examples

Two example problems are given:
* `SimpleOFE3`: A patch test problem using three meshes.
* `AMOREBracket`: A plane stress bracket with a hole.

Warning: The plane sweep algorithm is sensitive to rounding errors. For large problems, give coodinates small perturbations before computing the mesh overlay.

## Dependencies

* [SciPy](https://www.scipy.org)
* [NumPy](https://numpy.org)
* [matplotlab](https://matplotlib.org)
* [MATLAB](https://www.mathworks.com/products/matlab.html) (for post-processing only)

## Reference
J Huang, KJ Bathe. Overlapping finite element meshes in AMORE. <em>Advances in Engineering Software</em> 2020, 144:102791. [[journal]](https://doi.org/10.1016/j.advengsoft.2020.102791)

## License
[MIT](https://choosealicense.com/licenses/mit/)