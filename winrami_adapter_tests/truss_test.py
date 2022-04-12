from typing import List

from winrami_adapter.control import WinramiControl
from winrami_adapter_tests.dummy_truss import DummyTruss


STUB = 300


def test_create_truss():

    truss = DummyTruss()

    control = WinramiControl()

    top_chord_left_member_index = create_left_side_top_chord(control, truss)
    assert top_chord_left_member_index == 0

    top_chord_right_member_index = create_right_side_top_chord(control, truss)
    assert top_chord_right_member_index == 1

    bottom_chord_member_index = create_bottom_chord(control, truss)
    assert bottom_chord_member_index == 2

    bracing_member_indices = create_bracing_members(
        control=control,
        truss=truss,
        top_chord_left_member_index=top_chord_left_member_index,
        top_chord_right_member_index=top_chord_right_member_index,
        bottom_chord_member_index=bottom_chord_member_index
    )

    for index in range(6):
        assert bracing_member_indices[index * 2] == 3 + index * 2
        assert bracing_member_indices[index * 2 + 1] == 4 + index * 2

    # TODO Iiro: ristikon painon vertailu Winramiin.

    # TODO Ilari: control._native_model -> pitäisi tässä kohtaa olla valmis ja näyttää jotakuinkin järkevältä


def create_left_side_top_chord(control: WinramiControl, truss: DummyTruss) -> int:
    """
    Adds the left side top chord to calculation model and defines buckling length
    """

    top_chord_left_start_gp = truss.get_top_chord_left_gp()
    top_chord_left_start_gp_index = control.add_geometry_point_with_coordinates(
        top_chord_left_start_gp.x,
        top_chord_left_start_gp.y
    )
    top_chord_left_start_gp.add_wr_index(top_chord_left_start_gp_index)

    top_chord_left_end_gp = truss.get_ridge_gp()
    top_chord_left_end_gp_index = control.add_geometry_point_with_coordinates(
        top_chord_left_end_gp.x,
        top_chord_left_end_gp.y
    )
    top_chord_left_end_gp.add_wr_index(top_chord_left_end_gp_index)

    top_chord_left_member_index = control.add_member(
        top_chord_left_start_gp_index,
        top_chord_left_end_gp_index
    )

    control.unselect_all()
    control.select_member(top_chord_left_member_index)

    member = truss.get_top_chord_left_side_buckling()

    control.set_buckling_factor(abs(member.buckling_bend) * 1000, True)

    if abs(member.buckling_out) < 0.001: # ei nurjahdusta kohtisuorassa
        control.set_buckling_factor(1, False)
    elif member.buckling_out > 0.0:  # määrätty pituus 
        control.set_out_of_plane_buckling_length(member.buckling_out)

    return top_chord_left_member_index

def create_right_side_top_chord(control: WinramiControl, truss: DummyTruss) -> int:
    """
    Adds the right side top chord to calculation model and defines buckling length
    """

    top_chord_right_start_gp = truss.get_top_chord_right_gp()
    top_chord_right_start_gp_index = control.add_geometry_point_with_coordinates(
        top_chord_right_start_gp.x,
        top_chord_right_start_gp.y
    )
    top_chord_right_start_gp.add_wr_index(top_chord_right_start_gp_index)

    top_chord_right_end_gp = truss.get_ridge_gp()
    top_chord_right_end_gp_index = control.add_geometry_point_with_coordinates(
        top_chord_right_end_gp.x,
        top_chord_right_end_gp.y
    )
    top_chord_right_end_gp.add_wr_index(top_chord_right_end_gp_index)

    top_chord_right_member_index = control.add_member(
        top_chord_right_start_gp_index,
        top_chord_right_end_gp_index
    )

    control.unselect_all()
    control.select_member(top_chord_right_member_index)

    member = truss.get_top_chord_right_side_buckling()

    control.set_buckling_factor(abs(member.buckling_bend) * 1000, True)

    if abs(member.buckling_out) < 0.001: # ei nurjahdusta kohtisuorassa
        control.set_buckling_factor(1, False)
    elif member.buckling_out > 0.0:  # määrätty pituus 
        control.set_out_of_plane_buckling_length(member.buckling_out)

    return top_chord_right_member_index


def create_bottom_chord(control: WinramiControl, truss: DummyTruss) -> int:
    """
    Adds the bottom chord to calculation model and defines buckling length
    """

    bottom_chord_start_gp = truss.get_bottom_chord_left_gp()
    bottom_chord_start_gp_index = control.add_geometry_point_with_coordinates(
        bottom_chord_start_gp.x - STUB,
        bottom_chord_start_gp.y
    )
    bottom_chord_start_gp.add_wr_index(bottom_chord_start_gp_index)

    bottom_chord_end_gp = truss.get_bottom_chord_right_gp()
    bottom_chord_end_gp_index = control.add_geometry_point_with_coordinates(
        bottom_chord_end_gp.x + STUB,
        bottom_chord_end_gp.y
    )
    bottom_chord_end_gp.add_wr_index(bottom_chord_end_gp_index)

    bottom_chord_member_index = control.add_member(
        bottom_chord_start_gp_index,
        bottom_chord_end_gp_index,
    )

    control.unselect_all()
    control.select_member(bottom_chord_member_index)

    member = truss.get_bottom_chord_buckling()

    control.set_buckling_factor(abs(member.buckling_bend) * 1000, True)

    if member.buckling_out < -0.001:  # kerroin 
        control.set_buckling_factor(-abs(member.buckling_out) * 1000, True)
    elif abs(member.buckling_out) < 0.001:  # ei nurjahdusta kohtisuorassa
        control.set_buckling_factor(1, False)
    elif member.buckling_out > 0.0:  # määrätty pituus 
        control.set_out_of_plane_buckling_length(member.buckling_out)

    return bottom_chord_member_index


def create_bracing_members(
    *, control: WinramiControl,
    truss: DummyTruss,
    top_chord_left_member_index: int,
    top_chord_right_member_index: int,
    bottom_chord_member_index: int
) -> List[int]:
    """
    Adds the brancing members to calculation model and defines buckling lengths
    """

    start_offset = 150

    bracing_left_start_gp = truss.get_next_bracing_left_start_gp()
    bracing_left_end_gp = truss.get_next_bracing_left_end_gp()
    bracing_right_start_gp = truss.get_next_bracing_right_start_gp()
    bracing_right_end_gp = truss.get_next_bracing_right_end_gp()

    bracing_member_indices: List[int] = []

    def create_first_bracing_members() -> None:
        if not (bracing_left_start_gp and bracing_left_end_gp and bracing_right_start_gp and bracing_right_end_gp):
            return

        first_bracing_left_start_gp_index = control.add_geometry_point_to_member(
            top_chord_left_member_index,
            start_offset
        )
        bracing_left_start_gp.add_wr_index(first_bracing_left_start_gp_index)

        first_bracing_left_end_gp_index = control.add_geometry_point_to_member(
            bottom_chord_member_index,
            truss.calculate_length_along_bottom_chord(bracing_left_end_gp) + STUB
        )
        bracing_left_end_gp.add_wr_index(first_bracing_left_end_gp_index)

        first_bracing_left_member_index = control.add_member(
            first_bracing_left_start_gp_index,
            first_bracing_left_end_gp_index
        )
        bracing_member_indices.append(first_bracing_left_member_index)

        first_bracing_right_start_gp_index = control.add_geometry_point_to_member(
            top_chord_right_member_index,
            start_offset
        )
        bracing_right_start_gp.add_wr_index(first_bracing_right_start_gp_index)

        first_bracing_right_end_gp_index = control.add_geometry_point_to_member(
            bottom_chord_member_index,
            truss.calculate_length_along_bottom_chord(bracing_right_end_gp) + STUB
        )
        bracing_right_end_gp.add_wr_index(first_bracing_right_end_gp_index)

        first_bracing_right_member_index = control.add_member(
            first_bracing_right_start_gp_index,
            first_bracing_right_end_gp_index
        )
        bracing_member_indices.append(first_bracing_right_member_index)

    def create_remaining_bracing_members() -> None:
        nonlocal bracing_left_start_gp, bracing_left_end_gp, bracing_right_start_gp, bracing_right_end_gp

        while True:
            bracing_left_start_gp = truss.get_next_bracing_left_start_gp(bracing_left_start_gp)
            bracing_left_end_gp = truss.get_next_bracing_left_end_gp(bracing_left_end_gp)
            bracing_right_start_gp = truss.get_next_bracing_right_start_gp(bracing_right_start_gp)
            bracing_right_end_gp = truss.get_next_bracing_right_end_gp(bracing_right_end_gp)

            if not (bracing_left_start_gp and bracing_left_end_gp and bracing_right_start_gp and bracing_right_end_gp):
                break

            if not bracing_left_start_gp.has_wr_indices:
                bracing_left_end_gp_index = control.add_geometry_point_to_member(
                    top_chord_left_member_index,
                    truss.calculate_length_along_top_chord_left_side(bracing_left_start_gp)
                )
            else:
                bracing_left_end_gp_index = control.add_duplicate_geometry_point(
                    bracing_left_start_gp.last_wr_index
                )

            bracing_left_start_gp.add_wr_index(bracing_left_end_gp_index)

            if not bracing_left_end_gp.has_wr_indices:
                bracing_left_start_gp_index = control.add_geometry_point_to_member(
                    bottom_chord_member_index,
                    truss.calculate_length_along_bottom_chord(bracing_left_end_gp) + STUB
                )
            else:
                bracing_left_start_gp_index = control.add_duplicate_geometry_point(
                    bracing_left_end_gp.last_wr_index
                )

            bracing_left_end_gp.add_wr_index(bracing_left_start_gp_index)

            bracing_left_member_index = control.add_member(
                bracing_left_start_gp_index,
                bracing_left_end_gp_index
            )
            bracing_member_indices.append(bracing_left_member_index)

            if not bracing_right_start_gp.has_wr_indices:
                bracing_right_end_gp_index = control.add_geometry_point_to_member(
                    top_chord_right_member_index,
                    truss.calculate_length_along_top_chord_right_side(bracing_right_start_gp)
                )
            else:
                bracing_right_end_gp_index = control.add_duplicate_geometry_point(
                    bracing_right_start_gp.last_wr_index
                )

            bracing_right_start_gp.add_wr_index(bracing_right_end_gp_index)

            if not bracing_right_end_gp.has_wr_indices:
                bracing_right_start_gp_index = control.add_geometry_point_to_member(
                    bottom_chord_member_index,
                    truss.calculate_length_along_bottom_chord(bracing_right_end_gp) + STUB
                )
            else:
                bracing_right_start_gp_index = control.add_duplicate_geometry_point(
                    bracing_right_end_gp.last_wr_index
                )

            bracing_right_end_gp.add_wr_index(bracing_right_start_gp_index)

            bracing_right_member_index = control.add_member(
                bracing_right_start_gp_index,
                bracing_right_end_gp_index
            )
            bracing_member_indices.append(bracing_right_member_index)

        control.unselect_all()

        for member_index in bracing_member_indices:
            control.select_member(member_index)

        bracing_member = truss.get_bracing_member_buckling()

        control.set_buckling_factor(abs(bracing_member.buckling_bend) * 1000, True)
        control.set_buckling_factor(-abs(bracing_member.buckling_out) * 1000, True)

    create_first_bracing_members()
    create_remaining_bracing_members()

    return bracing_member_indices
