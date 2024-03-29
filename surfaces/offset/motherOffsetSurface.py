from surfaces.offset.offsetSurface import OffsetSurface
from surfaces.motherSurface import MotherSurface


class MotherOffsetSurface(OffsetSurface):
    def __init__(self, motherSurface: MotherSurface, offset: float) -> None:
        super().__init__(motherSurface, offset)
