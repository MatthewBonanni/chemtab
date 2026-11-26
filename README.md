# chemtab

A chemistry table generator for CFD. chemtab computes the thermochemical
states of a combustion model (steady flamelets for flamelet/progress-variable
(FPV) tables first, then FGM, heat loss, pressure and more), integrates
presumed PDFs, places the table's nodes to meet an interpolation tolerance, and
writes the table in [Mallard](https://github.com/MatthewBonanni/mallard)'s
format.

**Status: design.** Nothing is implemented yet. The agreed design is in
[docs/design.md](docs/design.md).

## Planned use

```bash
pip install chemtab
chemtab init h2-air        # writes chemtab.toml for a sample case
chemtab run chemtab.toml   # solves the flamelets and writes the table
```

Python users get the same thing from `chemtab.build(config)`, which returns the
table as NumPy arrays with named axes.

## Acknowledgments

chemtab is inspired by [FPVgen](https://github.com/IhmeGroup/FPVgen) (Ihme
Group, Stanford University), which generates FPV tables from Cantera
counterflow flames.

## License

Apache License 2.0; see [LICENSE](LICENSE).
