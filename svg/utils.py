from svgpathtools import svg2paths

import numpy as np
import math

POLYGON_MAX_CG_DISPLACEMENT = 1
JOIN_TAB_SAMPLE_DISTANCE = 0.1


class PolygonNotCenteredException(Exception):
    pass


def getPolygonFromSVGFile(filename, sampleDistance):
    paths, *_ = svg2paths(filename)

    if len(paths) != 1:
        raise NotImplementedError

    path = paths[0]

    segments = list(filter(lambda s: s.length() > 0, path))
    segmentsLenght = list([s.length() for s in segments])

    samples = int(math.ceil(path.length()/sampleDistance))

    def getK(sample, iSegment):
        segmentDistance = sampleDistance*sample - \
            sum(segmentsLenght[0:iSegment])
        return segmentDistance / segments[iSegment].length()

    def getPoint(sample, iSegment):
        point = segments[iSegment].point(getK(sample, iSegment))
        return np.array([point.real, point.imag])

    iSegment = 0
    for sample in range(samples):
        while getK(sample, iSegment) > 1:
            iSegment += 1

        yield getPoint(sample, iSegment)

    yield getPoint(0, 0)


def getTabJoinPolygonFromSVGFile(filename, width):
    tabJoinUnion = np.array(list(getPolygonFromSVGFile(
        filename, JOIN_TAB_SAMPLE_DISTANCE)))

    cg_x = sum(tabJoinUnion[:, 0])/len(tabJoinUnion)
    cg_y = sum(tabJoinUnion[:, 1])/len(tabJoinUnion)

    if ((cg_x**2 + cg_y**2) > POLYGON_MAX_CG_DISPLACEMENT**2):
        raise PolygonNotCenteredException()

    scale = width / (max(tabJoinUnion[0]) - min(tabJoinUnion[0]))

    return np.multiply(tabJoinUnion, scale)