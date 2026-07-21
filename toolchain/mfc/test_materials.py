"""Unit tests for the external material-file mechanism (mfc/materials.py): schema
validation, the provenance release_status gate, and the Cantera-style search order.
All coefficients are synthetic."""

import copy

import pytest
import yaml

from mfc.common import MFCException
from mfc.materials import find_material, load_material, material_fluid

SYNTHETIC_JWL = {
    "material": {
        "name": "synthetic_products",
        "eos": {
            "family": "jwl",
            "parameters": {"A": 5.0e11, "B": 8.0e9, "R1": 4.5, "R2": 1.2, "omega": 0.3, "rho0": 1600.0},
        },
        "provenance": {"citation": "synthetic coefficients for regression testing", "release_status": "synthetic"},
    }
}


def write_material(dirpath, filename, data):
    filepath = dirpath / filename
    filepath.write_text(yaml.safe_dump(data), encoding="utf-8")
    return str(filepath)


def test_load_and_fluid_mapping(tmp_path):
    filepath = write_material(tmp_path, "products.yaml", SYNTHETIC_JWL)

    material = load_material(filepath)
    assert material["name"] == "synthetic_products"
    assert material["eos"]["family"] == "jwl"

    entries = material_fluid(2, filepath)
    assert entries["fluid_pp(2)%eos"] == "jwl"
    assert entries["fluid_pp(2)%jwl_A"] == 5.0e11
    assert entries["fluid_pp(2)%jwl_rho0"] == 1600.0
    assert len(entries) == 7


def test_yaml11_unsigned_exponents(tmp_path):
    """PyYAML parses 5.0e11 (no exponent sign) as a string; the loader must still read
    it as the number a YAML 1.2 author intended."""
    filepath = tmp_path / "products.yaml"
    filepath.write_text(
        "material:\n"
        "  name: unsigned_exponents\n"
        "  eos:\n"
        "    family: jwl\n"
        "    parameters: {A: 5.0e11, B: 8.0e9, R1: 4.5, R2: 1.2, omega: 0.3, rho0: 1600.0}\n"
        "  provenance:\n"
        "    citation: synthetic\n"
        "    release_status: synthetic\n",
        encoding="utf-8",
    )
    entries = material_fluid(1, str(filepath))
    assert entries["fluid_pp(1)%jwl_A"] == 5.0e11
    assert entries["fluid_pp(1)%jwl_B"] == 8.0e9


def test_optional_q_validated_but_not_mapped(tmp_path):
    data = copy.deepcopy(SYNTHETIC_JWL)
    data["material"]["eos"]["parameters"]["Q"] = 4.0e6
    filepath = write_material(tmp_path, "products.yaml", data)

    assert load_material(filepath)["eos"]["parameters"]["Q"] == 4.0e6
    assert "fluid_pp(1)%jwl_Q" not in material_fluid(1, filepath)

    data["material"]["eos"]["parameters"]["Q"] = "large"
    filepath = write_material(tmp_path, "bad_q.yaml", data)
    with pytest.raises(MFCException, match="must be a number"):
        load_material(filepath)


def test_search_order(tmp_path, monkeypatch):
    case_dir = tmp_path / "case"
    public_dir = tmp_path / "public"
    case_dir.mkdir()
    public_dir.mkdir()
    monkeypatch.setenv("MFC_MATERIALS_DIR", str(public_dir))

    case_data = copy.deepcopy(SYNTHETIC_JWL)
    case_data["material"]["name"] = "from_case_dir"
    write_material(case_dir, "products.yaml", case_data)
    public_data = copy.deepcopy(SYNTHETIC_JWL)
    public_data["material"]["name"] = "from_public_dir"
    write_material(public_dir, "products.yaml", public_data)

    # the case directory is searched before the configured public directory
    assert load_material("products.yaml", str(case_dir))["name"] == "from_case_dir"
    # without a case-directory hit, the public directory resolves the bare name
    assert load_material("products.yaml", str(tmp_path))["name"] == "from_public_dir"

    monkeypatch.delenv("MFC_MATERIALS_DIR")
    with pytest.raises(MFCException, match="not found"):
        find_material("products.yaml", str(tmp_path))


def test_schema_rejections(tmp_path):
    rejections = [
        (["material"], None, "top-level 'material'"),
        (["material", "name"], "", "non-empty string"),
        (["material", "eos", "family"], "tillotson", "unknown EOS family"),
        (["material", "eos", "parameters", "R1"], None, "requires parameter 'R1'"),
        (["material", "eos", "parameters", "A"], "big", "must be a number"),
        (["material", "provenance"], None, "'material: provenance'"),
        (["material", "provenance", "citation"], "", "non-empty string"),
    ]
    for path, value, match in rejections:
        data = copy.deepcopy(SYNTHETIC_JWL)
        node = data
        for key in path[:-1]:
            node = node[key]
        if value is None:
            del node[path[-1]]
        else:
            node[path[-1]] = value
        filepath = write_material(tmp_path, "reject.yaml", data)
        with pytest.raises(MFCException, match=match):
            load_material(filepath)

    data = copy.deepcopy(SYNTHETIC_JWL)
    data["material"]["eos"]["parameters"]["detonation_velocity"] = 7000.0
    filepath = write_material(tmp_path, "reject.yaml", data)
    with pytest.raises(MFCException, match="unknown parameter"):
        load_material(filepath)


def test_release_status_gate(tmp_path):
    for status in (None, "restricted", "export_controlled", ""):
        data = copy.deepcopy(SYNTHETIC_JWL)
        if status is None:
            del data["material"]["provenance"]["release_status"]
        else:
            data["material"]["provenance"]["release_status"] = status
        filepath = write_material(tmp_path, "gate.yaml", data)
        with pytest.raises(MFCException, match="not releasable"):
            load_material(filepath)

    for status in ("synthetic", "open_literature", "public"):
        data = copy.deepcopy(SYNTHETIC_JWL)
        data["material"]["provenance"]["release_status"] = status
        filepath = write_material(tmp_path, "gate.yaml", data)
        assert load_material(filepath)["provenance"]["release_status"] == status
