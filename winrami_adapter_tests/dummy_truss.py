from math import sqrt
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class GeometryPoint:

    index: int
    x: float
    y: float
    wr_indices: List[int] = field(default_factory=list)

    def add_wr_index(self, index: int) -> None:
        self.wr_indices.append(index)

    @property
    def has_wr_indices(self) -> bool:
        return len(self.wr_indices) > 0

    @property
    def last_wr_index(self) -> int:
        return self.wr_indices[-1]


@dataclass
class MemberBuckling:
    buckling_bend: float
    buckling_out: float


class DummyTruss:

    gps = {
        1: GeometryPoint(index=1, x=0.0, y=0.0),
        3: GeometryPoint(index=3, x=2085.801, y=81.920),
        5: GeometryPoint(index=5, x=4057.200, y=159.347),
        7: GeometryPoint(index=7, x=6028.600, y=236.773),

        9: GeometryPoint(index=9, x=8000.0, y=314.200),

        11: GeometryPoint(index=11, x=9971.400, y=236.773),
        13: GeometryPoint(index=13, x=11942.800, y=159.347),
        15: GeometryPoint(index=15, x=13914.199, y=81.920),
        16: GeometryPoint(index=16, x=16000.0, y=0.0),

        2: GeometryPoint(index=2, x=2085.801, y=-2085.8),
        4: GeometryPoint(index=4, x=4057.200, y=-2085.800),
        6: GeometryPoint(index=6, x=6028.600, y=-2085.800),
        8: GeometryPoint(index=8, x=8000.000, y=-2085.800),
        10: GeometryPoint(index=10, x=9971.400, y=-2085.800),
        12: GeometryPoint(index=12, x=11942.800, y=-2085.800),
        14: GeometryPoint(index=14, x=13914.199, y=-2085.8),
    }

    bracing_member_start_points_left = [1, 3, 5, 5, 7, 9]
    bracing_member_end_points_left = [2, 2, 2, 6, 6, 6]
    bracing_member_start_points_right = [16, 15, 13, 13, 11, 9]
    bracing_member_end_points_right = [14, 14, 14, 10, 10, 10]

    top_chord_left_member = MemberBuckling(buckling_bend=-0.90, buckling_out=-1500.0)
    top_chord_right_member = MemberBuckling(buckling_bend=-0.90, buckling_out=-1500.0)
    bottom_chord_member = MemberBuckling(buckling_bend=-0.90, buckling_out=-0.90)
    bracing_member = MemberBuckling(buckling_bend=-0.75, buckling_out=-0.75)

    def __init__(self) -> None:
        self._bracing_member_start_points_left_index = 0
        self._bracing_member_end_points_left_index = 0
        self._bracing_member_start_points_right_index = 0
        self._bracing_member_end_points_right_index = 0

    def get_top_chord_left_gp(self) -> GeometryPoint:
        return self.gps[1]

    def get_ridge_gp(self) -> GeometryPoint:
        return self.gps[9]

    def get_top_chord_right_gp(self) -> GeometryPoint:
        return self.gps[16]

    def get_bottom_chord_left_gp(self) -> GeometryPoint:
        return self.gps[2]

    def get_bottom_chord_right_gp(self) -> GeometryPoint:
        return self.gps[14]

    def get_next_bracing_left_start_gp(self, previous_gp: Optional[GeometryPoint]=None) -> Optional[GeometryPoint]:
        if not previous_gp:
            self._bracing_member_start_points_left_index = 0
            return self.gps[self.bracing_member_start_points_left[self._bracing_member_start_points_left_index]]
        else:
            self._bracing_member_start_points_left_index += 1

            if self._bracing_member_start_points_left_index >= len(self.bracing_member_start_points_left):
                self._bracing_member_start_points_left_index = 0
                return None

            return  self.gps[self.bracing_member_start_points_left[self._bracing_member_start_points_left_index]]

    def get_next_bracing_left_end_gp(self, previous_gp: Optional[GeometryPoint]=None) -> Optional[GeometryPoint]:
        if not previous_gp:
            self._bracing_member_end_points_left_index = 0
            return self.gps[self.bracing_member_end_points_left[self._bracing_member_end_points_left_index]]
        else:
            self._bracing_member_end_points_left_index + 1

            if self._bracing_member_end_points_left_index >= len(self.bracing_member_end_points_left):
                self._bracing_member_end_points_left_index = 0
                return None

            return self.gps[self.bracing_member_end_points_left[self._bracing_member_end_points_left_index]]

    def get_next_bracing_right_start_gp(self, previous_gp: Optional[GeometryPoint]=None) -> Optional[GeometryPoint]:
        if not previous_gp:
            self._bracing_member_start_points_right_index = 0
            return self.gps[self.bracing_member_start_points_right[self._bracing_member_start_points_right_index]]
        else:
            self._bracing_member_start_points_right_index += 1

            if self._bracing_member_start_points_right_index >= len(self.bracing_member_start_points_right):
                self._bracing_member_start_points_right_index = 0
                return None

            return self.gps[self.bracing_member_start_points_right[self._bracing_member_start_points_right_index]]

    def get_next_bracing_right_end_gp(self, previous_gp: Optional[GeometryPoint]=None) -> Optional[GeometryPoint]:
        if not previous_gp:
            self._bracing_member_end_points_right_index = 0
            return self.gps[self.bracing_member_end_points_right[self._bracing_member_end_points_right_index]]
        else:
            self._bracing_member_end_points_right_index += 1

            if self._bracing_member_end_points_right_index >= len(self.bracing_member_end_points_right):
                self._bracing_member_end_points_right_index = 0
                return None

            return self.gps[self.bracing_member_end_points_right[self._bracing_member_end_points_right_index]]

    def get_top_chord_left_side_buckling(self) -> MemberBuckling:
        return self.top_chord_left_member

    def get_top_chord_right_side_buckling(self) -> MemberBuckling:
        return self.top_chord_right_member

    def get_bottom_chord_buckling(self) -> MemberBuckling:
        return self.bottom_chord_member

    def get_bracing_member_buckling(self) -> MemberBuckling:
        return self.bracing_member

    def calculate_length_along_top_chord_left_side(self, gp: GeometryPoint) -> float:
        top_chord_left_start = self.get_top_chord_left_gp()
        return sqrt(
            pow(gp.x - top_chord_left_start.x, 2) + 
            pow(gp.y - top_chord_left_start.y, 2)
        )

    def calculate_length_along_top_chord_right_side(self, gp: GeometryPoint) -> float:
        top_chord_right_start = self.get_ridge_gp()
        return sqrt(
            pow(gp.x - top_chord_right_start.x, 2) + 
            pow(gp.y - top_chord_right_start.y, 2)
        )

    def calculate_length_along_bottom_chord(self, gp: GeometryPoint) -> float:
        bottom_chord_start = self.get_bottom_chord_left_gp()
        return gp.x - bottom_chord_start.x
