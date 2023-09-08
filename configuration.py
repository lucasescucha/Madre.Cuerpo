# Default unit: milimeters (use Xm to use meters instead of mm)
configurationData = dict(
    surface=dict(
        dimensions=dict(
            # X axis
            width=500,
            # Y axis
            depth=400,
            # Z axis
            height=400,
        ),
        shape=dict(
            zeroHeight=200,
        )
    ),
    manufacture=dict(
        grid=dict(
            # X axis in arc lenght units (max value)
            width=5,
            # Y axis in arc lenght units (max value)
            height=5,
        ),
        mold=dict(
            puzzle=dict(
                thickness=4,
                pieces=dict(
                    # X axis in arc lenght units (max value)
                    width=150,
                    # Y axis in arc lenght units (max value)
                    height=150,
                ),
                tabs=dict(
                    # By default tab is y-axis oriented
                    file="svg/union.svg",
                    width=25,
                    baseThickness=2,
                    flangeSize=1
                ),
            ),
            panels=dict(
                dimensions=dict(
                    # unit: pieces
                    width=3,   
                    # unit: pieces
                    height=2,
                ),
                walls=dict(
                    height=30,
                    thickness=5,
                    unionMarks=3,
                    unionMarksDiameter=5,
                    screws=3,
                    screwsDiameter=3.1,
                    screwsHeadDiameter=5.7,
                    insertNutMargin=3,
                    insertNutDiameter=4,
                    insertNutDepth=4.5,
                    flangeSize=15,
                    flangeScrewsDiameter=3.1,
                    flangeScrewsHeadDiameter=5.7,
                    flangeNutSDimension=5.5
                ),
            ),
        ),
        steelFrame=dict(
            steel=dict(
                thickness=0.002,
                heigh=0.02,
            ),
            maxDistanceError=0.005,
            notchWidth=0.0006,
            notchLengt=0.019,
            simulateNotchAngle=False,
            divisions=2,
        ),
        fiberGlass=dict(
            firstLayerThickness=0.005,
        ),
    )
)