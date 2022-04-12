from winrami_adapter.control import BucklingMode, WinramiControl


def test_should_add_new_gp():
    control = WinramiControl()

    assert control.add_geometry_point_with_coordinates(0, 0) == 0
    assert control.add_geometry_point_with_coordinates(100, 100) == 1


def test_should_add_new_member():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)

    assert control.add_member(start_point, end_point) == 0


def test_should_select_added_elements():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    control.add_member(start_point, end_point)

    assert len(control.selected_points) == 2
    assert len(control.selected_members) == 1
    assert control.selected_points == {0, 1,}
    assert control.selected_members == {0}


def test_should_unselect_all_items():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    control.add_member(start_point, end_point)

    control.unselect_all()

    assert len(control.selected_points) == 0
    assert len(control.selected_members) == 0


def test_should_select_member():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    member = control.add_member(start_point, end_point)

    control.unselect_all()
    control.select_member(member)

    assert len(control.selected_points) == 0
    assert len(control.selected_members) == 1
    assert control.selected_members == {0}


def test_should_add_geometry_point_along_existing_member():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    member = control.add_member(start_point, end_point)

    assert control.add_geometry_point_to_member(member, 50) == 2


def test_should_add_duplicate_geometry_point():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    member = control.add_member(start_point, end_point)

    middle_point = control.add_geometry_point_to_member(member, 50)

    assert control.add_duplicate_geometry_point(middle_point) == 3


def test_should_update_buckling_mode():
    control = WinramiControl()

    start_point = control.add_geometry_point_with_coordinates(0, 0)
    end_point = control.add_geometry_point_with_coordinates(100, 100)
    member = control.add_member(start_point, end_point)

    control.unselect_all()
    control.select_member(member)

    assert control.current_buckling_mode == BucklingMode.IN_PLANE

    control.set_buckling_factor(0.9, True)

    assert control.current_buckling_mode == BucklingMode.OUT_OF_PLANE

    control.set_buckling_factor(1, False)

    assert control.current_buckling_mode == BucklingMode.IN_PLANE

    start_point = control.add_geometry_point_to_member(member, 50)
    end_point = control.add_geometry_point_with_coordinates(100, 0)
    member = control.add_member(start_point, end_point)

    control.unselect_all()
    control.select_member(member)

    assert control.current_buckling_mode == BucklingMode.IN_PLANE

    control.set_buckling_factor(0.9, True)

    assert control.current_buckling_mode == BucklingMode.OUT_OF_PLANE

    control.unselect_all()

    assert control.current_buckling_mode == BucklingMode.IN_PLANE