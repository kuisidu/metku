from dataclasses import dataclass
from enum import Enum, auto
from math import cos, sin, sqrt, tanh
from typing import Any, Dict, Set


@dataclass
class GP:

    x: float
    y: float


@dataclass
class Member:

    """Reference to starting geometry point"""
    start_gp: GP
    """Reference to ending geometry point"""
    end_gp: GP
    """Reference to native member ojbect"""
    native_member: Any # TODO Ilari: 'Any' tyypin voi muuttaa esim. 'SteelFrameMember' (tai mikä onkaan oikea luokka...)
                       #      niin ohjelman tarjoa intelisense on parempi

    @property
    def length(self) -> float:
        return sqrt(
            pow(self.end_gp.x - self.start_gp.x, 2) +
            pow(self.end_gp.y - self.start_gp.y, 2)
        )


class BucklingMode(Enum):

    IN_PLANE = auto()
    OUT_OF_PLANE = auto()


class WinramiControl:

    def __init__(self) -> None:
        self._next_gp_index = 0
        self._next_member_index = 0
        self._buckling_mode = BucklingMode.IN_PLANE

        self._geometry_points: Dict[int, GP] = {}
        self._members: Dict[int, Member] = {}
        self._selected_points: Set[int] = set()
        self._selected_members: Set[int] = set()

        # TODO Ilari: lisää Raamin laskentamalli:
        # self._native_model = Raami() (tai mikä onkaan oikea luokka)

    def add_geometry_point_with_coordinates(self, x: float, y: float) -> int:
        """
        Add a geometry point to specified coordinates.

        Returns an unique index for the new geometry point.
        """

        index = self._get_new_gp_index()
        self._geometry_points[index] = GP(x=x, y=y)
        self._selected_points.add(index)

        return index

    def add_geometry_point_to_member(self, member_index: int, offset_from_start: float) -> int:
        """
        Add a geometry point that is located along another member with offset
        length that is specified from the starting point of the member.

        Returns an unique index for the new geometry point.
        """

        member = self._members[member_index]

        if member.end_gp.x > member.start_gp.x:
            member_angle = tanh(
                (member.end_gp.y - member.start_gp.y) / (member.end_gp.x - member.start_gp.x)
            )

            x = cos(member_angle) * offset_from_start
            y = sin(member_angle) * offset_from_start
        else:
            member_angle = tanh(
                (member.end_gp.y - member.start_gp.y) / abs(member.end_gp.x - member.start_gp.x)
            )

            x = -cos(member_angle) * offset_from_start
            y = sin(member_angle) * offset_from_start

        return self.add_geometry_point_with_coordinates(x, y)

    def add_duplicate_geometry_point(self, index_of_previous_point: int) -> int:
        """
        Add a geometry point that is located in the same coordinates as an existing geometry point.

        Returns an unique index for the new geometry point.
        """

        index = self._get_new_gp_index()

        self._geometry_points[index] = GP(
            x=self._geometry_points[index_of_previous_point].x,
            y=self._geometry_points[index_of_previous_point].y
        )
        self._selected_points.add(index)

        return index

    def add_member(self, start_gp_index: int, end_gp_index: int) -> int:
        """
        Add member between two geometry points. Member spans from start point to end point.

        Returns an unique index for the new member.
        """

        start_point = self._geometry_points[start_gp_index]
        end_point = self._geometry_points[end_gp_index]

        native_member = None  # TODO Ilari: Luo Raamin käyttämä osa jollain oletus profiililla ja nurjahduksilla.
        #                            Osan alkupisteen x ja y löytyy: start_point.x ja start_point.y.
        #                            Osan loppupisteen x ja y löytyy: end_point.x ja end_point.y.
        #                            Huom! koordinaatit voivat olla negatiivisia. Onko ongelma?

        index = self._get_new_member_index()

        self._members[index] = Member(
            start_gp=start_point,
            end_gp=end_point,
            native_member=native_member
        )
        self._selected_members.add(index)

        return index

    def unselect_all(self) -> None:
        """
        Unselect all items.
        """

        self._selected_points.clear()
        self._selected_members.clear()
        self._buckling_mode = BucklingMode.IN_PLANE

    def select_member(self, member_index: int) -> None:
        """
        Select member with an index.
        """

        self._selected_members.add(member_index)

    def set_buckling_factor(self, buckling_factor: float, buckling_enabled: bool) -> None:
        """
        On the first call:
            - If buckling_enabled is True: set the in-plane buckling factor for selected members.
            - If buckling_enabled is False: disabled in-plane buckling for selected members (should not happen).

        On the second call:
            - If buckling_enabled is True: set the out-of-plane buckling factor for selected members.
            - If buckling_enabled is False: disabled out-of-plane buckling for selected members.
        """

        factor = buckling_factor if buckling_enabled else None

        for index in self._selected_members:
            member = self._members[index]

            if self._buckling_mode == BucklingMode.IN_PLANE:
                # TODO Ilari:
                #   - Mikäli factor ei ole None -> aseta vahvan(?) suunnan nurjahduskerroin
                #     member.native_member.XXX
                #   - Mikäli factor on None -> aseta eie nurjahdusta vahvassa(?) sunnassa
                #     (tätä tilannetta ei pitäisi FrameCalc optimoinnissa ikinä tulla, onko edes mielekäs?):
                #     member.native_member.XXX
                pass
            elif self._buckling_mode == BucklingMode.OUT_OF_PLANE:
                # TODO Ilari:
                #   - Mikäli factor ei ole None -> aseta heikon(?) suunnan nurjahduskerroin
                #     member.native_member....
                #   - Mikäli factor on None -> aseta ei nurjahdusta heikossa(?) sunnassa:
                #     (tämä tilanne voi tulla FrameCalc optimoinnissa):
                #     member.native_member.XXX
                pass

        self._buckling_mode = BucklingMode.IN_PLANE \
            if self._buckling_mode == BucklingMode.OUT_OF_PLANE \
            else BucklingMode.OUT_OF_PLANE

    def set_out_of_plane_buckling_length(self, buckling_length: float) -> None:
        """
        Set the out-of-plane buckling length for selected members.
        """

        for index in self._selected_members:
            member = self._members[index]
            factor = buckling_length / member.length

            # TODO Ilari: aseta heikon(?) suunnan nurjahduskerroin:
            #      member.native_member.XXX
            #      Ps. onko kertoimen laskenta noin ihan ok?

    def calculate(self) -> None:
        """Calculate model."""

        # TODO Ilari: tässä kohtaa lasketaan (ainakin) FEM ja ratkaistaan kuormitukset:
        #      self._native_model.structural_analysis(load_id, "REM")
        #      Tässä kohtaa halutaan siis laskea kaikki kuormitustapaukset; vaatiko erityisiä toimenpiteitä?
        #      Täytyy myöhemmin miettiä tarviiko tässä kohtaa laskea myös mitoitus vai
        #      voidaanko se tehdä erikseen vain tarpeen mukaan.

    @property
    def selected_points(self) -> Set[int]:
        """Return copy of selected point indices. Used for testing."""
        return self._selected_points.copy()

    @property
    def selected_members(self) -> Set[int]:
        """Return copy of selected member indices. Used for testing."""
        return self._selected_members.copy()

    @property
    def current_buckling_mode(self) -> BucklingMode:
        """Return currently active buckling mode. Used for testing."""
        return self._buckling_mode

    def _get_new_gp_index(self) -> int:
        index = self._next_gp_index
        self._next_gp_index += 1
        return index

    def _get_new_member_index(self) -> int:
        index = self._next_member_index
        self._next_member_index += 1
        return index
