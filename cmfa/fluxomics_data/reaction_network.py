from operator import gt, lt
from typing import List

from pydantic import BaseModel, Field, computed_field, model_validator


class Compound(BaseModel):
    id: str
    n_labellable_atoms: int = Field(ge=1)


class Transition(BaseModel):
    compound_id: str
    atom_pattern: str
    coef: float


class Reaction(BaseModel):
    id: str
    transitions: List[Transition]

    @model_validator(mode="after")
    def check_atom_balance(self) -> "Reaction":
        atoms_in_str, atoms_out_str = (
            "".join(
                sorted(t.atom_pattern for t in self.transitions if f(t.coef, 0))
            )
            for f in (lt, gt)
        )
        msg = (
            f"Reaction {self.id} has unbalanced atom transitions."
            f"\n\tAtoms in: {atoms_in_str}"
            f"\n\tAtoms out: {atoms_out_str}\n"
        )
        assert atoms_in_str == atoms_out_str, msg
        return self


class ReactionNetwork(BaseModel):
    id: str
    reactions: List[Reaction]

    @computed_field
    def compounds(self) -> List[Compound]:
        out = []
        for r in self.reactions:
            for t in r.transitions:
                new = Compound(
                    id=t.compound_id, n_labellable_atoms=len(t.atom_pattern)
                )
                for cpd in out:
                    if cpd.id == new.id:
                        assert cpd.n_labellable_atoms == new.n_labellable_atoms
                if not any(t.compound_id == c.id for c in out):
                    out += [new]
        return out

    @model_validator(mode="after")
    def check_compounds(self) -> "ReactionNetwork":
        return self
