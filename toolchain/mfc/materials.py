"""External, provenance-tagged material files.

Material files follow the Cantera mechanism-file model: the repository carries the
equation family, while calibrated coefficient sets live outside it and are referenced
by path or name. A material file is YAML with one material:

    material:
      name: user_material
      eos:
        family: jwl
        parameters:
          A: 5.0e11
          B: 8.0e9
          R1: 4.5
          R2: 1.2
          omega: 0.3
          rho0: 1600.0
          Q: 0.0            # optional; reserved for reactive models, unused here
      provenance:
        citation: where the numbers come from
        release_status: synthetic

Search locations mirror Cantera mechanisms: the spec as an explicit path, then the
case directory, then the configured public material directory ($MFC_MATERIALS_DIR).
Analytic EOS parameters remain runtime values, so loading a material never triggers a
rebuild; tabular EOS data stays external by the same rule.
"""

import os

import yaml

from .common import MFCException

# Per-family required parameter keys; optional keys are accepted and validated but
# carried only in the material dict, not mapped into fluid_pp entries.
EOS_FAMILIES = {
    "jwl": {"required": ("A", "B", "R1", "R2", "omega", "rho0"), "optional": ("Q",)},
}

# release_status values that may enter a case; anything else is refused so that
# restricted calibrations cannot flow into the solver by accident.
RELEASABLE_STATUSES = ("synthetic", "open_literature", "public")


def find_material(spec: str, case_dirpath: str = None) -> str:
    candidates = [spec]
    if case_dirpath is not None:
        candidates.append(os.path.join(case_dirpath, spec))
    if os.environ.get("MFC_MATERIALS_DIR"):
        candidates.append(os.path.join(os.environ["MFC_MATERIALS_DIR"], spec))

    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate

    raise MFCException(f"Material file '{spec}' not found. Searched: {', '.join(candidates)}.")


def _require(condition: bool, filepath: str, message: str) -> None:
    if not condition:
        raise MFCException(f"Material file '{filepath}': {message}")


def load_material(spec: str, case_dirpath: str = None) -> dict:
    """Locate, parse, and validate a material file; returns the inner material dict."""
    filepath = find_material(spec, case_dirpath)
    with open(filepath, "rt", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    _require(isinstance(data, dict) and isinstance(data.get("material"), dict), filepath, "top-level 'material' mapping is missing.")
    material = data["material"]

    name = material.get("name")
    _require(isinstance(name, str) and name.strip() != "", filepath, "'material: name' must be a non-empty string.")

    eos = material.get("eos")
    _require(isinstance(eos, dict), filepath, "'material: eos' mapping is missing.")
    family = eos.get("family")
    _require(family in EOS_FAMILIES, filepath, f"unknown EOS family '{family}'. Supported: {', '.join(EOS_FAMILIES)}.")

    parameters = eos.get("parameters")
    _require(isinstance(parameters, dict), filepath, "'eos: parameters' mapping is missing.")
    required, optional = EOS_FAMILIES[family]["required"], EOS_FAMILIES[family]["optional"]
    for key in required:
        _require(key in parameters, filepath, f"family '{family}' requires parameter '{key}'.")
    for key in list(parameters):
        _require(key in required or key in optional, filepath, f"unknown parameter '{key}' for family '{family}'.")
        if isinstance(parameters[key], str):
            # PyYAML reads YAML 1.1, whose float exponents need a sign; accept 5.0e11 too
            try:
                parameters[key] = float(parameters[key])
            except ValueError:
                pass
        value = parameters[key]
        _require(isinstance(value, (int, float)) and not isinstance(value, bool), filepath, f"parameter '{key}' must be a number.")

    provenance = material.get("provenance")
    _require(isinstance(provenance, dict), filepath, "'material: provenance' mapping is missing.")
    citation = provenance.get("citation")
    _require(isinstance(citation, str) and citation.strip() != "", filepath, "'provenance: citation' must be a non-empty string.")
    release_status = provenance.get("release_status")
    _require(
        release_status in RELEASABLE_STATUSES,
        filepath,
        f"release_status '{release_status}' is not releasable into a case. Allowed: {', '.join(RELEASABLE_STATUSES)}.",
    )

    return material


def material_fluid(index: int, spec: str, case_dirpath: str = None) -> dict:
    """Return the fluid_pp(index)%... case-dictionary entries for a material file.

    Only the EOS selection and its analytic parameters are mapped; reactive fields
    such as Q are validated but not consumed by nonreactive runs.
    """
    material = load_material(spec, case_dirpath)
    family = material["eos"]["family"]
    parameters = material["eos"]["parameters"]
    prefix = f"fluid_pp({index})%"

    entries = {prefix + "eos": family}
    for key in EOS_FAMILIES[family]["required"]:
        entries[f"{prefix}{family}_{key}"] = float(parameters[key])
    return entries
